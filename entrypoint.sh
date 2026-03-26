#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# SMART — Somatic Mutation Annotation and Reporting Tool
# ==============================================================================
# Dockerised entrypoint — replaces the original main.sh, removing all
# apptainer/singularity calls since tools are installed in the container.
# ==============================================================================

usage() {
cat <<EOF
SMART — Somatic Mutation Annotation and Reporting Tool (Dockerised)

Usage:
  docker run -v /path/to/data:/data -v /path/to/refs:/refs smart \\
    <ONCOKB_TOKEN> --transcripts-file FILE --ref-dir /refs [OPTIONS]

Required:
  <ONCOKB_TOKEN>           OncoKB API token
  --transcripts-file FILE  Transcript list for prioritisation
  --ref-dir DIR            Reference resources directory (mounted at /refs)

  The ref-dir should contain:
    DIR/liftover/hg19ToHg38.over.chain
    DIR/liftover/hg38.fa
    DIR/Plugins
    DIR/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz
    DIR/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz
    DIR/REVEL/new_tabbed_revel_grch38.tsv.gz
    DIR/ClinVar/clinvar.vcf.gz
    DIR/CIVIC/civic_01_10_25.vcf.gz
    DIR/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz
    DIR/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz

Options:
  --pass / --no-pass               Enable/disable PASS filtering (default: ON)
  --liftover / --no-liftover       Enable/disable liftover (default: ON)
  --clean-tmp / --keep-tmp         Delete/keep intermediate files (default: clean)
  --clean-tables / --keep-tables   Delete/keep table files after post_analysis (default: clean)
  --help                           Show this help

Volumes:
  Mount your data directory (containing OriginalVcf/) at /data
  Mount your reference directory at /refs (or wherever --ref-dir points)

Examples:
  docker run --rm \\
    -v \$(pwd):/data \\
    -v /mnt/data1/vep_test:/refs \\
    smart \\
    "\$ONCOKB_TOKEN" \\
    --transcripts-file /data/TSO500_transcripts_list.txt \\
    --ref-dir /refs

EOF
}

if [[ $# -lt 1 ]] || [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    usage
    exit 0
fi

ONCOKB_TOKEN="$1"
shift

# Defaults
DO_PASS_FILTER=1
DO_LIFTOVER=1
CLEAN_TMP=1
CLEAN_TABLES=1
TRANSCRIPTS_FILE=""
REF_DIR=""

# Parse flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-pass)        DO_PASS_FILTER=0; shift ;;
        --pass)           DO_PASS_FILTER=1; shift ;;
        --no-liftover)    DO_LIFTOVER=0;    shift ;;
        --liftover)       DO_LIFTOVER=1;    shift ;;
        --clean-tmp)      CLEAN_TMP=1;      shift ;;
        --keep-tmp)       CLEAN_TMP=0;      shift ;;
        --clean-tables)   CLEAN_TABLES=1;   shift ;;
        --keep-tables)    CLEAN_TABLES=0;   shift ;;
        --transcripts-file)
            [[ $# -lt 2 ]] && { echo "ERROR: --transcripts-file requires a path"; exit 2; }
            TRANSCRIPTS_FILE="$2"; shift 2 ;;
        --ref-dir)
            [[ $# -lt 2 ]] && { echo "ERROR: --ref-dir requires a path"; exit 2; }
            REF_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; usage; exit 2 ;;
    esac
done

# Validate required arguments
[[ -z "$TRANSCRIPTS_FILE" ]] && { echo "ERROR: --transcripts-file is required"; usage; exit 2; }
[[ ! -f "$TRANSCRIPTS_FILE" ]] && { echo "ERROR: Transcript file not found: $TRANSCRIPTS_FILE"; exit 2; }
[[ -z "$REF_DIR" ]] && { echo "ERROR: --ref-dir is required"; usage; exit 2; }
[[ ! -d "$REF_DIR" ]] && { echo "ERROR: Reference directory not found: $REF_DIR"; exit 2; }

# Reference files
CHAIN="$REF_DIR/liftover/hg19ToHg38.over.chain"
REF="$REF_DIR/liftover/hg38.fa"
VEP_DIR="$REF_DIR"
VEP_PLUGIN_DIR="$REF_DIR/Plugins"
SPLICEAI_SNV="$REF_DIR/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$REF_DIR/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
REVEL_FILE="$REF_DIR/REVEL/new_tabbed_revel_grch38.tsv.gz"
CLINVAR_VCF="$REF_DIR/ClinVar/clinvar.vcf.gz"
CIVIC_VCF="$REF_DIR/CIVIC/civic_01_10_25.vcf.gz"
LOEUF_FILE="$REF_DIR/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
CANCER_HOTSPOTS_VCF="$REF_DIR/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"

# Validate reference files
REQUIRED_PATHS=(
    "$CHAIN" "$REF" "$VEP_PLUGIN_DIR"
    "$SPLICEAI_SNV" "$SPLICEAI_INDEL" "$REVEL_FILE"
    "$CLINVAR_VCF" "$CIVIC_VCF" "$LOEUF_FILE" "$CANCER_HOTSPOTS_VCF"
)
for path in "${REQUIRED_PATHS[@]}"; do
    [[ ! -e "$path" ]] && { echo "ERROR: Required resource not found: $path"; exit 2; }
done

echo "============================================================"
echo "SMART — Somatic Mutation Annotation and Reporting Tool"
echo "============================================================"
echo "Tool versions:"
echo "  VEP:              $(vep --help 2>&1 | grep -oP 'ensembl-vep\s*:\s*\K\S+' || echo 'unknown')"
echo "  GATK:             $(gatk --version 2>&1 | grep -oP '\d+\.\d+\.\d+\.\d+' | head -1 || echo 'unknown')"
echo "  bcftools:          $(bcftools --version | head -1 || echo 'unknown')"
echo "  Python:            $(python3 --version 2>&1 || echo 'unknown')"
echo ""
echo "Config:"
echo "  PASS filter:      $([[ $DO_PASS_FILTER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Liftover:         $([[ $DO_LIFTOVER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Clean tmp:        $([[ $CLEAN_TMP -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Clean tables:     $([[ $CLEAN_TABLES -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Transcript file:  $TRANSCRIPTS_FILE"
echo "  Reference dir:    $REF_DIR"
echo "============================================================"
echo

SCRIPT_DIR="/opt/smart/scripts"

# Create output directories
INPUT_DIR="./OriginalVcf"
FILTERED_DIR="./FilteredVcf"
OUTPUT_DIR="./LiftOverVcf"
REJECT_DIR="./LiftOverVcf/Rejected"
ANNOTATED_DIR="./AnnotatedVcf"
ONCOKB_DIR="./OncoKB_VCF"
TABLE_DIR="./Table"
FINAL_TABLE_DIR="./FINAL_Table"

mkdir -p "$FILTERED_DIR" "$OUTPUT_DIR" "$REJECT_DIR" "$ANNOTATED_DIR" "$ONCOKB_DIR" "$TABLE_DIR" "$FINAL_TABLE_DIR"

echo "Cleaning previous results..."
rm -f "$OUTPUT_DIR"/*gz "$OUTPUT_DIR"/*tbi "$OUTPUT_DIR"/*.vcf 2>/dev/null || true
rm -f "$REJECT_DIR"/* 2>/dev/null || true
rm -f "$ANNOTATED_DIR"/* 2>/dev/null || true
rm -f "$FILTERED_DIR"/* 2>/dev/null || true
rm -f "$ONCOKB_DIR"/* 2>/dev/null || true
rm -f "$TABLE_DIR"/* 2>/dev/null || true
rm -f "$FINAL_TABLE_DIR"/* 2>/dev/null || true
echo "Cleanup complete."

COUNT_FILE="./variant_counts.txt"
echo -e "Sample\tOriginal_Variants\tPASS_Variants\tLifted_Variants\tVEP_Annotated_Variants\tOncoKB_Annotated_Variants\tVariantsInTable" > "$COUNT_FILE"

for vcf in "$INPUT_DIR"/*.vcf.gz; do
    [[ -e "$vcf" ]] || continue

    basename=$(basename "$vcf")
    sample="${basename%.dragen.concat_snv_sv.vcf.gz}"

    echo "###################### Processing: $sample ##########################################"

    # --- PASS filtering ---
    FILTERED_VCF="$FILTERED_DIR/${sample}.PASS_only.vcf.gz"
    FILTERED_VCF_TBI="${FILTERED_VCF}.tbi"

    if [[ $DO_PASS_FILTER -eq 1 ]]; then
        echo ">>> PASS filtering: $sample"
        zcat "$vcf" \
            | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 == "PASS"' \
            | bgzip > "$FILTERED_VCF"
        tabix -p vcf "$FILTERED_VCF"
        INPUT_FOR_NEXT="$FILTERED_VCF"
    else
        echo ">>> PASS filtering DISABLED: $sample"
        INPUT_FOR_NEXT="$vcf"
    fi

    # --- LiftOver ---
    LIFTED_VCF_GZ="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf.gz"
    LIFTED_VCF="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf"
    LIFTED_VCF_TBI="${LIFTED_VCF_GZ}.tbi"
    REJECTED_VCF="$REJECT_DIR/${sample}.hg19tohg38_rejected.vcf.gz"
    REJECTED_VCF_TBI="${REJECTED_VCF}.tbi"

    if [[ $DO_LIFTOVER -eq 1 ]]; then
        echo ">>> LiftOver: $sample"
        # Direct GATK call (no apptainer)
        gatk LiftoverVcf \
            -C "$CHAIN" \
            -I "$INPUT_FOR_NEXT" \
            -O "$LIFTED_VCF_GZ" \
            -R "$REF" \
            --WRITE_ORIGINAL_POSITION \
            --REJECT "$REJECTED_VCF"
        gunzip -c "$LIFTED_VCF_GZ" > "$LIFTED_VCF"
        VEP_INPUT="$LIFTED_VCF"
        echo ">>> LiftOver complete: $sample"
    else
        echo ">>> LiftOver DISABLED: $sample"
        VEP_INPUT="$INPUT_FOR_NEXT"
    fi

    # --- VEP Annotation (direct call, no apptainer) ---
    echo ">>> VEP annotation: $sample"
    ANNOTATED_VCF="$ANNOTATED_DIR/${sample}_annotated.vcf"
    vep \
        -i "$VEP_INPUT" \
        -o "$ANNOTATED_VCF" \
        --vcf --everything --offline \
        --dir "$VEP_DIR" \
        --dir_plugins "$VEP_PLUGIN_DIR" \
        --assembly GRCh38 \
        --no_escape \
        --fasta "$REF" \
        --plugin SpliceAI,snv="$SPLICEAI_SNV",indel="$SPLICEAI_INDEL",cutoff=0.5 \
        --plugin REVEL,"$REVEL_FILE" \
        --custom "$CLINVAR_VCF",ClinVar,vcf,exact,0,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNSIGSCV,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ONCDN,ONCDNINCL,ONCDISDB,ONCDISDBINCL,ONC,ONCINCL,ONCREVSTAT,ONCSCV,ONCCONF,ORIGIN,RS,SCIDN,SCIDNINCL,SCIDISDB,SCIDISDBINCL,SCIREVSTAT,SCI,SCIINCL,SCISCV \
        --custom "$CIVIC_VCF",CIViC,vcf,exact,0,GN,VT,CSQ \
        --plugin LOEUF,file="$LOEUF_FILE",match_by=transcript \
        --custom "$CANCER_HOTSPOTS_VCF",CancerHotspots,vcf,exact,0,HOTSPOT,HOTSPOT_GENE,HOTSPOT_HGVSp,HOTSPOT3D,HOTSPOT3D_GENE,HOTSPOT3D_HGVSp,HOTSPOTNC,HOTSPOTNC_GENE,HOTSPOTNC_HGVSc

    # --- OncoKB annotation ---
    echo ">>> OncoKB annotation: $sample"
    ONCOKB_VCF="$ONCOKB_DIR/${sample}.oncoKB.vcf"
    python3 "$SCRIPT_DIR/oncokb2.0.py" "$ANNOTATED_VCF" "$ONCOKB_VCF" \
        --oncokb_token "$ONCOKB_TOKEN" \
        --tumor_mode generic \
        --preferred-transcripts "$TRANSCRIPTS_FILE"

    # --- VCF to Table ---
    echo ">>> VCF2TABLE: $sample"
    TABLE_CSV="$TABLE_DIR/${sample}.oncoKB.csv"
    python3 "$SCRIPT_DIR/vcf2table.py" \
        --vcf "$ONCOKB_VCF" \
        --transcripts "$TRANSCRIPTS_FILE" \
        --output "$TABLE_CSV" \
        --debug

    # --- MafAnnotator ---
    echo ">>> MafAnnotator: $sample"
    ANNOTATED_MAF="$FINAL_TABLE_DIR/${sample}.oncoKB.maf"
    python3 /opt/oncokb-annotator/MafAnnotator.py \
        -i "$TABLE_CSV" \
        -o "$ANNOTATED_MAF" \
        -b "$ONCOKB_TOKEN"

    # --- Variant counts ---
    ORIGINAL_COUNT=$(zgrep -vc "^#" "$vcf")
    if [[ $DO_PASS_FILTER -eq 1 ]]; then
        PASS_COUNT=$(zgrep -vc "^#" "$FILTERED_VCF")
    else
        PASS_COUNT="$ORIGINAL_COUNT"
    fi
    if [[ $DO_LIFTOVER -eq 1 ]]; then
        LIFTED_COUNT=$(grep -vc "^#" "$LIFTED_VCF")
    else
        LIFTED_COUNT="$PASS_COUNT"
    fi
    ANNOTATED_COUNT=$(grep -vc "^#" "$ANNOTATED_VCF")
    ONCOKB_COUNT=$(grep -vc "^#" "$ONCOKB_VCF")
    TABLE_COUNT=$(( $(wc -l < "$TABLE_CSV") - 1 ))

    echo -e "${sample}\t${ORIGINAL_COUNT}\t${PASS_COUNT}\t${LIFTED_COUNT}\t${ANNOTATED_COUNT}\t${ONCOKB_COUNT}\t${TABLE_COUNT}" >> "$COUNT_FILE"

    # --- Cleanup ---
    if [[ $CLEAN_TMP -eq 1 ]]; then
        echo ">>> Cleaning temporary files for: $sample"
        rm -f "$FILTERED_VCF" "$FILTERED_VCF_TBI"
        rm -f "$LIFTED_VCF_GZ" "$LIFTED_VCF" "$LIFTED_VCF_TBI"
        rm -f "$REJECTED_VCF" "$REJECTED_VCF_TBI"
        rm -f "$ANNOTATED_VCF"
        rm -f "$ONCOKB_VCF"
    fi

    echo "--- $sample complete ---"
done

echo "###################### Post Analysis ##########################################"
python3 "$SCRIPT_DIR/post_analysis.py"

if [[ $CLEAN_TABLES -eq 1 ]]; then
    echo ">>> Cleaning table files after post_analysis"
    rm -f "$TABLE_DIR"/*
    rm -f "$FINAL_TABLE_DIR"/*
fi

echo "============================================================"
echo "SMART PIPELINE COMPLETE"
echo "============================================================"
