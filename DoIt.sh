echo "Cleaning previous results..."

rm -f -r ./LiftOverVcf/*
rm -f ./AnnotatedVcf/*
rm -f ./Stats/*
rm -f ./FilteredVcf/*
rm -f ./OncoKB_VCF/*
rm -f ./Table/*

BCFTOOLS="apptainer exec --bind $(pwd):/mnt --bind /mnt/data1:/mnt/data1 --bind /mnt/scratch1:/mnt/scratch1 docker://quay.io/biocontainers/bcftools:1.21--h3a4d415_1 bcftools"

echo "Cleanup complete. Generating new results..."

INPUT_DIR="./OriginalVcf"
FILTERED_DIR="./FilteredVcf"
OUTPUT_DIR="./LiftOverVcf"
REJECT_DIR="./LiftOverVcf/Rejected"
ANNOTATED_DIR="./AnnotatedVcf"
STATS="./Stats"
ONCOKB_DIR="./OncoKB_VCF"
TABLE="./Table"

CHAIN="/root/vep_sample_test/manuel/DoItAll/liftover/hg19ToHg38.over.chain"
REF="/root/vep_sample_test/manuel/DoItAll/liftover/hg38.fa"

mkdir -p "$FILTERED_DIR" "$OUTPUT_DIR" "$REJECT_DIR" "$ANNOTATED_DIR" "$STATS" "$ONCOKB_DIR" "$TABLE"

COUNT_FILE="./variant_counts.txt"
echo -e "Sample\tOriginal\tPASS\tLifted\tVEP_Annotated\tOncoKB_Annotated" > "$COUNT_FILE"


for vcf in "$INPUT_DIR"/*.dragen.concat_snv_sv.vcf.gz; do
    [ -e "$vcf" ] || continue

    basename=$(basename "$vcf")
    sample="${basename%.dragen.concat_snv_sv.vcf.gz}"

    echo "###################### Processing: $sample ##########################################"

    # Count total variants in original VCF
    ORIGINAL_COUNT=$(zgrep -v "^#" "$vcf" | wc -l)

    #$BCFTOOLS stats -f PASS $vcf > $STATS/${sample}_PASS.stat.txt


    ###############################################
    # Step 1 — PASS filtering
    ###############################################
    FILTERED_VCF="$FILTERED_DIR/${sample}.PASS_only.vcf.gz"

    zcat "$vcf" \
    | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 == "PASS"' \
    | bgzip > "$FILTERED_VCF"

    tabix -p vcf "$FILTERED_VCF"

    #PASS_COUNT=$(zgrep -v "^#" "$FILTERED_VCF" | wc -l)


    echo "###################### LiftOver STARTING: $sample ##########################################"
    ###############################################
    # Step 2 — LiftOver
    ###############################################
    LIFTED_VCF_GZ="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf.gz"
    LIFTED_VCF="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf"
    REJECTED_VCF="$REJECT_DIR/${sample}.hg19tohg38_rejected.vcf.gz"

    apptainer exec --bind "$(pwd)" --bind /mnt:/mnt \
        docker://broadinstitute/gatk:4.6.0.0 \
        /gatk/gatk LiftoverVcf \
        -C "$CHAIN" \
        -I "$FILTERED_VCF" \
        -O "$LIFTED_VCF_GZ" \
        -R "$REF" \
        --WRITE_ORIGINAL_POSITION \
        --REJECT "$REJECTED_VCF"

    gunzip -c $LIFTED_VCF_GZ > $LIFTED_VCF

    LIFTED_COUNT=$(grep -v "^#" "$LIFTED_VCF" | wc -l)
    echo "###################### LiftOver END: $sample ##########################################"

    echo "###################### VEP annotating STARTING: $sample ##########################################"
    ###############################################
    # Step 3 — VEP Annotation
    ###############################################
    ANNOTATED_VCF="$ANNOTATED_DIR/${sample}_annotated.vcf"

    /root/ensembl-vep/vep \
        -i "$LIFTED_VCF" \
        -o "$ANNOTATED_VCF" \
        --vcf --everything --offline \
        --dir /mnt/data1/vep_test \
        --dir_plugins /mnt/data1/vep_test/Plugins \
        --assembly GRCh38 \
        --fasta "$REF" \
        --plugin SpliceAI,snv=/mnt/data1/vep_test/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/mnt/data1/vep_test/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
        --plugin REVEL,/mnt/data1/vep_test/REVEL/new_tabbed_revel_grch38.tsv.gz \
        --custom /mnt/data1/vep_test/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNSIGSCV,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ONCDN,ONCDNINCL,ONCDISDB,ONCDISDBINCL,ONC,ONCINCL,ONCREVSTAT,ONCSCV,ONCCONF,ORIGIN,RS,SCIDN,SCIDNINCL,SCIDISDB,SCIDISDBINCL,SCIREVSTAT,SCI,SCIINCL,SCISCV \
        --custom /mnt/data1/vep_test/CIVIC/civic_01_10_25.vcf.gz,CIViC,vcf,exact,0,GN,VT,CSQ \
        --plugin LOEUF,file=/mnt/data1/vep_test/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz,match_by=transcript \
        --custom /mnt/data1/vep_test/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz,CancerHotspots,vcf,exact,0,HOTSPOT,HOTSPOT_GENE,HOTSPOT_HGVSp,HOTSPOT3D,HOTSPOT3D_GENE,HOTSPOT3D_HGVSp,HOTSPOTNC,HOTSPOTNC_GENE,HOTSPOTNC_HGVSc

    echo "###################### VEP annotating END: $sample ##########################################"
    echo "###################### OncoKB annotation STARTING: $sample ##########################################"
    ###############################################
    # Step 4 — OncoKB annotation
    ###############################################
    ONCOKB_VCF="$ONCOKB_DIR/${sample}.oncoKB.vcf"

    python3 oncokb2.0.py "$ANNOTATED_VCF" "$ONCOKB_VCF" --tumor_mode=generic
    echo "###################### OncoKB annotation END: $sample ##########################################"

    #ONCOKB_COUNT=$(grep -v "^#" "$ONCOKB_VCF" | wc -l)
    ###############################################
    # Step 4.5 — Convert Annotated VCF to Table
    ###############################################
    echo "###################### VCF2TABLE: $sample ##########################################"
    TABLE_CSV="$TABLE/${sample}.oncoKB.csv"

    python3 vcf2table.py \
        --vcf "$ONCOKB_VCF" \
        --transcripts TSO500_transcripts_list.txt \
        --output "$TABLE_CSV" \
        --debug

    ###############################################
    # Step 5 — Count summary
    ###############################################
    PASS_COUNT=$(zgrep -v "^#" "$FILTERED_VCF" | wc -l)
    ANNOTATED_COUNT=$(grep -v "^#" "$ANNOTATED_VCF" | wc -l)
    ONCOKB_COUNT=$(grep -v "^#" "$ONCOKB_VCF" | wc -l)

    echo -e "${sample}\t${ORIGINAL_COUNT}\t${PASS_COUNT}\t${LIFTED_COUNT}\t${ANNOTATED_COUNT}\t${ONCOKB_COUNT}" >> "$COUNT_FILE"

    echo "---"
done



# ANNOTATED_DIR="./AnnotatedVcf"
# ONCOKB_DIR="./OncoKB_VCF"
# TABLE="./Table"

# mkdir -p "$ONCOKB_DIR" "$TABLE"

# for vcf in "$ANNOTATED_DIR"/*_annotated.vcf; do
#     [ -e "$vcf" ] || continue

#     basename=$(basename "$vcf")
#     sample="${basename%_annotated.vcf}"

#     echo "###################### Processing: $sample ##########################################"

#     ###############################################
#     # Step 4 — OncoKB annotation
#     ###############################################
#     ONCOKB_VCF="$ONCOKB_DIR/${sample}.oncoKB.vcf"

#     echo "###################### OncoKB annotation STARTING: $sample ##########################################"
#     python3 oncokb2.0.py "$vcf" "$ONCOKB_VCF" --tumor_mode=generic
#     echo "###################### OncoKB annotation END: $sample ##########################################"

#     ###############################################
#     # Step 4.5 — Convert Annotated VCF to Table
#     ###############################################
#     TABLE_CSV="$TABLE/${sample}.oncoKB.csv"

#     echo "###################### VCF2TABLE: $sample ##########################################"
#     python3 vcf2table.py \
#         --vcf "$ONCOKB_VCF" \
#         --transcripts TSO500_transcripts_list.txt \
#         --output "$TABLE_CSV" \
#         --debug

#     echo "---"
# done