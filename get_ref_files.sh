#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Clinical Genomics Pipeline – Reference Data Installer
# ==============================================================================
#
# This script downloads and prepares the reference datasets required to run
# the variant annotation pipeline.
#
# Features:
#   - Idempotent: existing files are skipped
#   - Creates indexes where needed
#   - Creates a pipeline_config.sh file in the current directory
#   - Records paths and resource versions where possible
#
# Notes:
#   - Some resources (Google Drive / BaseSpace) may require manual validation
#   - This script assumes wget, git, gunzip, zgrep, realpath, apptainer are available
#   - tabix is run on the host system in this version
#
# ==============================================================================

echo "================================================================================"
echo " Clinical Genomics Pipeline – Reference Data Installer"
echo ""
echo "This script downloads and prepares all reference datasets required to run"
echo "the variant annotation pipeline."
echo ""
echo "Resources installed include:"
echo "  • Human reference genome (GRCh38)"
echo "  • Genome indexes required by GATK and samtools"
echo "  • hg19 → GRCh38 liftover chain file"
echo "  • Public annotation databases (ClinVar, CIViC, gnomAD constraints, etc.)"
echo "  • VEP cache and plugins for offline annotation"
echo ""
echo "The installer is idempotent:"
echo "  • Existing files are detected and skipped"
echo "  • Only missing resources are downloaded"
echo "  • The script can safely be re-run if interrupted"
echo ""
echo "All resources will be installed in the directory:"
echo "  ref_files/"
echo "================================================================================"

REF_ROOT="ref_files"
mkdir -p "$REF_ROOT"

# -------------------------
# Container wrappers
# -------------------------
echo "[INIT] Setting container wrappers"

SAMTOOLS="apptainer exec \
    --bind /mnt/data1:/mnt/data1 \
    --bind /mnt/scratch1:/mnt/scratch1 \
    --bind $(pwd):$(pwd) \
    docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 samtools"

GATK="apptainer exec \
    --bind /mnt:/mnt \
    --bind $(pwd):$(pwd) \
    docker://broadinstitute/gatk:4.1.3.0"

VEP="apptainer exec \
    --bind /mnt:/mnt \
    --bind $(pwd):$(pwd) \
    docker://ensemblorg/ensembl-vep"

# -------------------------
# Helper functions
# -------------------------
log_step () {
  echo ""
  echo "================================================================================"
  echo "$1"
  echo "================================================================================"
}

download_if_missing () {
  # $1 = URL, $2 = OUTPUT FILE
  local url="$1"
  local out="$2"

  if [[ -f "$out" ]]; then
    echo "[SKIP] Already exists: $out"
    return 0
  fi

  echo "[DL] $url"
  echo "     -> $out"

  if ! wget -O "$out" "$url"; then
    echo "[ERROR] Download failed: $url"
    rm -f "$out"
    exit 1
  fi
}

mkdir_if_missing () {
  local d="$1"
  if [[ -d "$d" ]]; then
    echo "[OK] Directory exists: $d"
  else
    echo "[MKDIR] Creating: $d"
    mkdir -p "$d"
  fi
}

# -------------------------
# STEP 1: Liftover chain
# -------------------------
log_step "[STEP 1] Liftover chain (hg19 → hg38)"

mkdir_if_missing "$REF_ROOT/liftover"

CHAIN_GZ="$REF_ROOT/liftover/hg19ToHg38.over.chain.gz"
CHAIN="$REF_ROOT/liftover/hg19ToHg38.over.chain"

if [[ -f "$CHAIN" ]]; then
  echo "[SKIP] Liftover chain already exists: $CHAIN"
else
  download_if_missing \
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    "$CHAIN_GZ"

  echo "[RUN] Decompressing chain file -> $CHAIN"
  gunzip -c "$CHAIN_GZ" > "$CHAIN"
fi

# -------------------------
# STEP 2: GRCh38 reference + indexes
# -------------------------
log_step "[STEP 2] GRCh38 reference genome (hg38.fa + .fai + .dict)"

mkdir_if_missing "$REF_ROOT/reference"

REF_FA_GZ="$REF_ROOT/reference/hg38.fa.gz"
REF_FA="$REF_ROOT/reference/hg38.fa"
REF_FAI="$REF_ROOT/reference/hg38.fa.fai"
REF_DICT="$REF_ROOT/reference/hg38.dict"

if [[ -f "$REF_FA" ]]; then
  echo "[SKIP] Reference FASTA already exists: $REF_FA"
else
  download_if_missing \
    "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz" \
    "$REF_FA_GZ"

  echo "[RUN] Decompressing reference FASTA -> $REF_FA"
  gunzip -f "$REF_FA_GZ"
fi

if [[ -f "$REF_FAI" ]]; then
  echo "[SKIP] FASTA index already exists: $REF_FAI"
else
  echo "[RUN] Creating FASTA index (.fai) with samtools"
  $SAMTOOLS faidx "$REF_FA"
fi

if [[ -f "$REF_DICT" ]]; then
  echo "[SKIP] Sequence dictionary already exists: $REF_DICT"
else
  echo "[RUN] Creating sequence dictionary (.dict) with GATK"
  $GATK /gatk/gatk CreateSequenceDictionary \
      -R "$REF_FA" \
      -O "$REF_DICT"
fi

# -------------------------
# STEP 3: ClinVar
# -------------------------
log_step "[STEP 3] ClinVar (GRCh38 VCF + tabix index)"

mkdir_if_missing "$REF_ROOT/ClinVar"
CLINVAR_VCF="$REF_ROOT/ClinVar/clinvar.vcf.gz"
CLINVAR_TBI="$REF_ROOT/ClinVar/clinvar.vcf.gz.tbi"

download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
  "$CLINVAR_VCF"

if [[ -f "$CLINVAR_TBI" ]]; then
  echo "[SKIP] ClinVar tabix index already exists: $CLINVAR_TBI"
else
  echo "[RUN] Indexing ClinVar with tabix"
  tabix -p vcf "$CLINVAR_VCF"
fi

# -------------------------
# STEP 4: CIViC
# -------------------------
log_step "[STEP 4] CIViC (nightly download)"

mkdir_if_missing "$REF_ROOT/CIVIC"
CIVIC_FILE="$REF_ROOT/CIVIC/nightly-VariantSummaries.tsv"

download_if_missing \
  "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv" \
  "$CIVIC_FILE"

# -------------------------
# STEP 5: REVEL
# -------------------------
log_step "[STEP 5] REVEL (zip download)"

mkdir_if_missing "$REF_ROOT/REVEL"
REVEL_ZIP="$REF_ROOT/REVEL/revel_all_chromosomes.zip"

download_if_missing \
  "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" \
  "$REVEL_ZIP"

# -------------------------
# STEP 6: SpliceAI
# -------------------------
log_step "[STEP 6] SpliceAI (SNV + INDEL VCFs)"

mkdir_if_missing "$REF_ROOT/SpliceAI"

SPLICEAI_SNV="$REF_ROOT/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$REF_ROOT/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"

download_if_missing "https://basespace.illumina.com/s/otSPW8hnhaZR" "$SPLICEAI_SNV"
download_if_missing "https://basespace.illumina.com/s/jI4aHqRKoGqE" "$SPLICEAI_INDEL"

log_step "WARNING: This will FAIL if you dont have an Illumina account. This is not available outside https://basespace.illumina.com/"

cp /mnt/data1/vep_test/SpliceAI/* "$REF_ROOT/SpliceAI"

# -------------------------
# STEP 7: gnomAD constraints (LOEUF)
# -------------------------
log_step "[STEP 7] gnomAD constraint metrics (LOEUF)"

mkdir_if_missing "$REF_ROOT/gnomAD_constraints"

LOEUF_GRCH38="$REF_ROOT/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_GRCH38_TBI="$REF_ROOT/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1pJMBwZQQLB2K4DelrWCEILzSadmLVybq" \
  "$LOEUF_GRCH38"

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1ZFg9d2VW_PC3VcsaSS84LaKeudw41PT8" \
  "$LOEUF_GRCH38_TBI"

# -------------------------
# STEP 8: VEP cache
# -------------------------
log_step "[STEP 8] VEP cache (homo_sapiens, GRCh38)"

VEP_CACHE_DIR="$REF_ROOT/homo_sapiens"

if [[ -d "$VEP_CACHE_DIR" ]] && [[ "$(ls -A "$VEP_CACHE_DIR" 2>/dev/null || true)" != "" ]]; then
  echo "[SKIP] VEP cache directory already present and non-empty: $VEP_CACHE_DIR"
else
  echo "[RUN] Installing VEP cache into: $REF_ROOT"
  $VEP /root/miniconda3/pkgs/ensembl-vep-88.9-0/bin/vep_install \
    -a cf \
    -s homo_sapiens \
    -y GRCh38 \
    -c "$REF_ROOT"
fi

# -------------------------
# STEP 9: VEP plugins
# -------------------------
log_step "[STEP 9] VEP Plugins"

PLUGINS_DIR="$REF_ROOT/Plugins"
if [[ -d "$PLUGINS_DIR/.git" ]]; then
  echo "[SKIP] VEP plugins repo already cloned: $PLUGINS_DIR"
else
  if [[ -d "$PLUGINS_DIR" ]]; then
    echo "[WARN] $PLUGINS_DIR exists but is not a git repo. Skipping clone to avoid overwriting."
    echo "       If you want to re-clone, remove the directory and rerun."
  else
    echo "[RUN] Cloning Ensembl VEP plugins repository..."
    git clone https://github.com/Ensembl/VEP_plugins.git "$PLUGINS_DIR"
  fi
fi

# -------------------------
# STEP 10: Cancer Hotspots
# -------------------------
log_step "[STEP 10] Download Cancer Hotspots (hg38)"

mkdir_if_missing "$REF_ROOT/CancerHotSpots"

HOTSPOT_VCF="$REF_ROOT/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"

if [[ -f "$HOTSPOT_VCF" ]]; then
  echo "[SKIP] Cancer Hotspots VCF already exists: $HOTSPOT_VCF"
else
  echo "[DL] Downloading Cancer Hotspots dataset..."
  if ! wget -O "$HOTSPOT_VCF" \
    "https://raw.githubusercontent.com/charlottekyng/cancer_hotspots/master/resources/hg38.hotspots_changv2_gao_nc.vcf.gz"; then
    echo "[ERROR] Download failed: Cancer Hotspots"
    rm -f "$HOTSPOT_VCF"
    exit 1
  fi
fi

# -------------------------
# STEP 11: Write config file
# -------------------------
log_step "[STEP 11] Writing pipeline configuration file"

CONFIG_FILE="pipeline_config.sh"

# Resource versions / metadata
INSTALL_DATE="$(date +%Y-%m-%d)"
REFERENCE_VERSION="UCSC_hg38"
CHAIN_VERSION="UCSC_hg19ToHg38"
CIVIC_VERSION="civic_01_10_25"
REVEL_VERSION="v1.3"
SPLICEAI_VERSION="SpliceAIv1.3" 
LOEUF_VERSION="custom_grch38"
CANCER_HOTSPOTS_VERSION="charlottekyng/cancer_hotspots master"
VEP_CACHE_VERSION="GRCh38"
VEP_PLUGINS_VERSION="git clone"

CLINVAR_VERSION="$(zgrep -m1 '^##fileDate=' "$CLINVAR_VCF" 2>/dev/null | cut -d= -f2 || true)"
if [[ -z "$CLINVAR_VERSION" ]]; then
  CLINVAR_VERSION="unknown"
fi

cat > "$CONFIG_FILE" <<EOF
# --------------------------------------------------
# Clinical Genomics Pipeline configuration
# Generated automatically by install_references.sh
# --------------------------------------------------

INSTALL_DATE="$INSTALL_DATE"
REF_ROOT="$(realpath "$REF_ROOT")"

REF_GENOME="\$REF_ROOT/reference/hg38.fa"
REF_GENOME_VERSION="$REFERENCE_VERSION"

CHAIN_FILE="\$REF_ROOT/liftover/hg19ToHg38.over.chain"
CHAIN_FILE_VERSION="$CHAIN_VERSION"

CLINVAR_VCF="\$REF_ROOT/ClinVar/clinvar.vcf.gz"
CLINVAR_VERSION="$CLINVAR_VERSION"

CIVIC_FILE="\$REF_ROOT/CIVIC/nightly-VariantSummaries.tsv"
CIVIC_VERSION="$CIVIC_VERSION"

REVEL_FILE="\$REF_ROOT/REVEL/revel_all_chromosomes.zip"
REVEL_VERSION="$REVEL_VERSION"

SPLICEAI_SNV="\$REF_ROOT/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="\$REF_ROOT/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
SPLICEAI_VERSION="$SPLICEAI_VERSION"

LOEUF_FILE="\$REF_ROOT/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_FILE_INDEX="\$REF_ROOT/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"
LOEUF_VERSION="$LOEUF_VERSION"

VEP_CACHE_DIR="\$REF_ROOT/homo_sapiens"
VEP_CACHE_VERSION="$VEP_CACHE_VERSION"

VEP_PLUGINS_DIR="\$REF_ROOT/Plugins"
VEP_PLUGINS_VERSION="$VEP_PLUGINS_VERSION"

CANCER_HOTSPOTS_VCF="\$REF_ROOT/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
CANCER_HOTSPOTS_VERSION="$CANCER_HOTSPOTS_VERSION"
EOF

echo "[OK] Configuration file created: $CONFIG_FILE"

echo ""
echo "================================================================================"
echo " DONE"
echo "================================================================================"
echo "Resources installed under: $(realpath "$REF_ROOT")"
echo "Configuration file written to: $(realpath "$CONFIG_FILE")"
echo "If something failed, rerun this script—completed steps will be skipped."







