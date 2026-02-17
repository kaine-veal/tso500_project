
set -euo pipefail

#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------
# This pipeline processes raw DRAGEN VCF files through a full clinical‑style annotation
# workflow. This script is designed to work across many samples (VCFs) 
# It performs the following steps for each sample:
#
#   1. Cleanup of previous results to ensure a fresh run.
#       Temporary files are not deleted until new execution. For debugging, etc
#   2.  Filter variants to retain only PASS calls.
#
#   3.  LiftOver variants from hg19 to GRCh38.
#
#   4.  Annotate lifted variants using Ensembl VEP with
#      multiple plugins (SpliceAI, REVEL, LOEUF, ClinVar, CIViC, CancerHotspots, etc.).
#
#   5. Annotate the VEP‑processed VCF using OncoKB (oncokb2.0.py), producing an
#      OncoKB‑enriched VCF for each sample.
#
#   6. Convert the OncoKB‑annotated VCF into a tabular format using vcf2table.py.
#      This step extracts VEP annotations, INFO fields, FORMAT/sample fields, and
#      expands the embedded OncoKB JSON into individual columns.
#
#   7. Run the OncoKB MafAnnotator to produce a final MAF‑style annotated table
#      suitable for downstream reporting or clinical interpretation.
#
#   8. Generate a summary file (variant_counts2.txt) containing per‑sample counts
#      for each major processing stage.
#
# How to run:
#   - Conda activate TSO500
#   - Put Original Dragen VCF files in AnnotatedVCF folder
#   - nohup ./DoIt.sh > DoIt.log 2>&1 &
#   - This will run in the wglsbi03 background and a log file will be created.
#   - or run directly ./DoIt.sh 
# Author: Manuel, Kaine and Ian
# -----------------------------------------------------------------------------------------

#ONCOKB_TOKEN="b1c06077-370f-47d9-b5df-7a9fede45da1"
ONCOKB_TOKEN="$1"

# Create the folders where tmp and final files will be saved
INPUT_DIR="./OriginalVcf"
FILTERED_DIR="./FilteredVcf"
OUTPUT_DIR="./LiftOverVcf"
REJECT_DIR="./LiftOverVcf/Rejected"
ANNOTATED_DIR="./AnnotatedVcf"
ONCOKB_DIR="./OncoKB_VCF"
TABLE="./Table"
FINAL_TABLE="./FINAL_Table"
mkdir -p "$FILTERED_DIR" "$OUTPUT_DIR" "$REJECT_DIR" "$ANNOTATED_DIR" "$ONCOKB_DIR" "$TABLE" "$FINAL_TABLE"


echo "Cleaning previous results..."
# Better to do this just in case some tools don't like rewrite

rm -f  ./LiftOverVcf/*gz ./LiftOverVcf/*tbi
rm -f  ./LiftOverVcf/Rejected/* 
rm -f ./AnnotatedVcf/*
rm -f ./FilteredVcf/*
rm -f ./OncoKB_VCF/*
rm -f ./Table/*
rm -f ./FINAL_Table/*

echo "Cleanup complete. Generating new results..."

# The chain and the fasta file needed for lift over
CHAIN="/root/vep_sample_test/manuel/DoItAll/liftover/hg19ToHg38.over.chain"
REF="/root/vep_sample_test/manuel/DoItAll/liftover/hg38.fa"

# When working with many samples, I wanted to ensure the number of variants
# in each step, to ensure I dont lose variants:
COUNT_FILE="./variant_counts.txt"

# The header of this table
echo -e "Sample\tOriginal_Variants\tPASS_Variants\tLifted_Variants\tVEP_Annotated_Variants\tOncoKB_Annotated_Variants\tVariantsInTable\tVariantsInFinalTable" > "$COUNT_FILE"

# The loop that does all the work.
# For each VCF found in $INPUT_DIR (the original VCF without annotation)
for vcf in "$INPUT_DIR"/*.dragen.concat_snv_sv.vcf.gz; do
    [ -e "$vcf" ] || continue

    basename=$(basename "$vcf")
    sample="${basename%.dragen.concat_snv_sv.vcf.gz}"

    echo "###################### Processing: $sample ##########################################"

    ###############################################
    # Step 2 — PASS filtering
    ###############################################

    # I filter variants to select only PASS variants
    # Becaouse low quality variants return many issues downstream
    # e.g., many Low QUAL variants are rejected by the lift over

    FILTERED_VCF="$FILTERED_DIR/${sample}.PASS_only.vcf.gz"

    zcat "$vcf" \
    | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 == "PASS"' \
    | bgzip > "$FILTERED_VCF"

    tabix -p vcf "$FILTERED_VCF"

    echo "###################### LiftOver STARTING: $sample ##########################################"
    ###############################################
    # Step 3 — LiftOver
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
        --no_escape \
        --fasta "$REF" \
        --plugin SpliceAI,snv=/mnt/data1/vep_test/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/mnt/data1/vep_test/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
        --plugin REVEL,/mnt/data1/vep_test/REVEL/new_tabbed_revel_grch38.tsv.gz \
        --custom /mnt/data1/vep_test/ClinVar/clinvar.vcf.gz,ClinVar,vcf,exact,0,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNSIGSCV,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ONCDN,ONCDNINCL,ONCDISDB,ONCDISDBINCL,ONC,ONCINCL,ONCREVSTAT,ONCSCV,ONCCONF,ORIGIN,RS,SCIDN,SCIDNINCL,SCIDISDB,SCIDISDBINCL,SCIREVSTAT,SCI,SCIINCL,SCISCV \
        --custom /mnt/data1/vep_test/CIVIC/civic_01_10_25.vcf.gz,CIViC,vcf,exact,0,GN,VT,CSQ \
        --plugin LOEUF,file=/mnt/data1/vep_test/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz,match_by=transcript \
        --custom /mnt/data1/vep_test/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz,CancerHotspots,vcf,exact,0,HOTSPOT,HOTSPOT_GENE,HOTSPOT_HGVSp,HOTSPOT3D,HOTSPOT3D_GENE,HOTSPOT3D_HGVSp,HOTSPOTNC,HOTSPOTNC_GENE,HOTSPOTNC_HGVSc

    # echo "###################### VEP annotating END: $sample ##########################################"
    echo "###################### OncoKB annotation STARTING: $sample ##########################################"
    ###############################################
    # Step 4 — OncoKB annotation
    ###############################################
    ONCOKB_VCF="$ONCOKB_DIR/${sample}.oncoKB.vcf"

    # Skip if OncoKB output already exists
    if [ -f "$ONCOKB_VCF" ]; then
        echo "OncoKB VCF already exists for $sample — skipping OncoKB step."
        continue
    fi

    python3 oncokb2.0.py "$ANNOTATED_VCF" "$ONCOKB_VCF" --tumor_mode=generic --oncokb_token=$ONCOKB_TOKEN
    # Generic mode does not assign a tumour‑specific annotation.
    # This is likely the recommended approach for this project, but if needed,
    # oncokb2.0.py can extract the tumour name from the filename.
    # However, the tumour name in the filename is not always valid, so some
    # replacements are applied to assign a standardised tumour name.

    echo "###################### OncoKB annotation END: $sample ##########################################"

    ###############################################
    # Step 4.5 — Convert Annotated VCF to Table
    ###############################################
    echo "###################### VCF2TABLE: $sample ##########################################"
    TABLE_CSV="$TABLE/${sample}.oncoKB.csv"
    # This script converts a VCF into a lightweight MAF-style file,
    # which is sufficient for compatibility with the OncoKB MafAnnotator.
    python3 vcf2table.py \
        --vcf "$ONCOKB_VCF" \
        --transcripts TSO500_transcripts_list.txt \
        --output "$TABLE_CSV" \
        --debug
    
    ###############################################
    # Step 4.6 — OncoKB MafAnnotator
    ###############################################
    echo "###################### MafAnnotator: $sample ##########################################"

    ANNOTATED_MAF="$FINAL_TABLE/${sample}.oncoKB.maf"

    python3 ../oncokb-annotator/MafAnnotator.py \
        -i "$TABLE_CSV" \
        -o "$ANNOTATED_MAF" \
        -b "$ONCOKB_TOKEN"


    # ###############################################
    # # Step 5 — Count summary
    # ###############################################
    
    ORIGINAL_COUNT=$(zgrep -v "^#" "$vcf" | wc -l)
    PASS_COUNT=$(zgrep -v "^#" "$FILTERED_VCF" | wc -l)
    LIFTED_COUNT=$(grep -v "^#" "$LIFTED_VCF" | wc -l)
    ANNOTATED_COUNT=$(grep -v "^#" "$ANNOTATED_VCF" | wc -l)
    ONCOKB_COUNT=$(grep -v "^#" "$ONCOKB_VCF" | wc -l)
    TABLE_COUNT=$(($(wc -l < "$TABLE_CSV") - 1))

    echo -e "${sample}\t${ORIGINAL_COUNT}\t${PASS_COUNT}\t${LIFTED_COUNT}\t${ANNOTATED_COUNT}\t${ONCOKB_COUNT}\t${TABLE_COUNT}" >> "$COUNT_FILE"

    echo "---"
done

echo "###################### Post Analysis step ##########################################"

# This is a final data management step to merge all tables in one
# and organise the final table deleting unnecesary columns, adding the
# the second header, etc
python3 post_analysis.py

echo "DONE"