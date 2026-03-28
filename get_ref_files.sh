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
#   - Verifies checksums where available
#   - Checks disk space before starting
#
# Notes:
#   - SpliceAI requires an Illumina BaseSpace account; files must be manually
#     staged to $SPLICEAI_STAGING_DIR before running this script
#   - This script assumes wget, git, gunzip, zgrep, realpath, apptainer,
#     and tabix are available on the host
#
# ==============================================================================

# ==============================================================================
# CONFIGURATION — edit these variables to match your environment
# ==============================================================================

# Base install directory (defaults to the running user's home directory)
BASE_DIR="${HOME}"

# Where reference files will be stored
REF_ROOT="${BASE_DIR}/ref_files"

# SpliceAI files must be manually downloaded from BaseSpace and placed here
# before running this script (BaseSpace requires an Illumina login)
SPLICEAI_STAGING_DIR="${BASE_DIR}/spliceai_staging"

# Minimum free disk space required to proceed (in GB)
MIN_DISK_GB=200

# ==============================================================================

echo "================================================================================"
echo " Clinical Genomics Pipeline – Reference Data Installer"
echo ""
echo " Running as user : $(whoami)"
echo " Home directory  : ${HOME}"
echo " Install root    : ${REF_ROOT}"
echo " SpliceAI source : ${SPLICEAI_STAGING_DIR}"
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
echo "================================================================================"

# -------------------------
# Preflight checks
# -------------------------

echo ""
echo "[PREFLIGHT] Checking available disk space on ${BASE_DIR}..."
AVAIL_GB=$(df -BG "${BASE_DIR}" | awk 'NR==2 {gsub("G",""); print $4}')
echo "[PREFLIGHT] Available: ${AVAIL_GB} GB  |  Required: ${MIN_DISK_GB} GB"

if (( AVAIL_GB < MIN_DISK_GB )); then
  echo "[ERROR] Insufficient disk space."
  echo "        Available : ${AVAIL_GB} GB"
  echo "        Required  : ${MIN_DISK_GB} GB"
  echo "        Free up space under ${BASE_DIR} and re-run."
  exit 1
fi

echo "[PREFLIGHT] Disk space OK."

# -------------------------
# Container wrappers
# -------------------------
echo ""
echo "[INIT] Setting container wrappers"

# In the Uni Server we have 
module load apptainer

# -------------------------
# Container wrappers
# -------------------------
echo "[INIT] Setting container wrappers"

SAMTOOLS_IMG="docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"

SAMTOOLS="apptainer exec --fakeroot \
    --bind ${BASE_DIR}:${BASE_DIR} \
    --bind $(pwd):$(pwd) \
    ${SAMTOOLS_IMG} samtools"

TABIX="apptainer exec --fakeroot \
    --bind ${BASE_DIR}:${BASE_DIR} \
    --bind $(pwd):$(pwd) \
    ${SAMTOOLS_IMG} tabix"

GATK="apptainer exec --fakeroot \
    --bind ${BASE_DIR}:${BASE_DIR} \
    --bind $(pwd):$(pwd) \
    docker://broadinstitute/gatk:4.1.3.0"

VEP="apptainer exec --fakeroot \
    --bind ${BASE_DIR}:${BASE_DIR} \
    --bind $(pwd):$(pwd) \
    docker://ensemblorg/ensembl-vep:release_114.0"


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
  # Usage: download_if_missing <URL> <OUTPUT_FILE> [EXPECTED_SHA256]
  local url="$1"
  local out="$2"
  local expected_sha="${3:-}"

  if [[ -f "$out" ]]; then
    echo "[SKIP] Already exists: $out"
    # Still verify checksum if one was provided
    if [[ -n "$expected_sha" ]]; then
      local actual_sha
      actual_sha=$(sha256sum "$out" | awk '{print $1}')
      if [[ "$actual_sha" != "$expected_sha" ]]; then
        echo "[ERROR] Checksum mismatch for existing file: $out"
        echo "        Expected : $expected_sha"
        echo "        Got      : $actual_sha"
        exit 1
      fi
      echo "[OK] Checksum verified: $out"
    fi
    return 0
  fi

  echo "[DL] $url"
  echo "     -> $out"

  local tmp_out="${out}.tmp"

  if ! wget --show-progress -O "$tmp_out" "$url"; then
    echo "[ERROR] Download failed: $url"
    rm -f "$tmp_out"
    exit 1
  fi

  # Verify checksum before moving into place
  if [[ -n "$expected_sha" ]]; then
    local actual_sha
    actual_sha=$(sha256sum "$tmp_out" | awk '{print $1}')
    if [[ "$actual_sha" != "$expected_sha" ]]; then
      echo "[ERROR] Checksum mismatch after download: $url"
      echo "        Expected : $expected_sha"
      echo "        Got      : $actual_sha"
      rm -f "$tmp_out"
      exit 1
    fi
    echo "[OK] Checksum verified."
  fi

  mv "$tmp_out" "$out"
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

mkdir_if_missing "${REF_ROOT}/liftover"

CHAIN_GZ="${REF_ROOT}/liftover/hg19ToHg38.over.chain.gz"
CHAIN="${REF_ROOT}/liftover/hg19ToHg38.over.chain"

if [[ -f "$CHAIN" ]]; then
  echo "[SKIP] Liftover chain already exists: $CHAIN"
else
  download_if_missing \
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    "$CHAIN_GZ"

  echo "[RUN] Decompressing chain file -> $CHAIN"
  gunzip -c "$CHAIN_GZ" > "$CHAIN"
  echo "[OK] Chain file ready: $CHAIN"
fi

# -------------------------
# STEP 2: GRCh38 reference + indexes
# -------------------------
log_step "[STEP 2] GRCh38 reference genome (hg38.fa + .fai + .dict)"

mkdir_if_missing "${REF_ROOT}/reference"

REF_FA_GZ="${REF_ROOT}/liftover/hg38.fa.gz"
REF_FA="${REF_ROOT}/liftover/hg38.fa"
REF_FAI="${REF_ROOT}/liftover/hg38.fa.fai"
REF_DICT="${REF_ROOT}/liftover/hg38.dict"

if [[ -f "$REF_FA" ]]; then
  echo "[SKIP] Reference FASTA already exists: $REF_FA"
else
  download_if_missing \
    "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz" \
    "$REF_FA_GZ"

  echo "[RUN] Decompressing reference FASTA -> $REF_FA"
  # Use -c (stdout) to preserve the .gz file, protecting against interruption
  gunzip -c "$REF_FA_GZ" > "$REF_FA"
  echo "[OK] Reference FASTA ready: $REF_FA"
fi

if [[ -f "$REF_FAI" ]]; then
  echo "[SKIP] FASTA index already exists: $REF_FAI"
else
  echo "[RUN] Creating FASTA index (.fai) with samtools"
  $SAMTOOLS faidx "$REF_FA"
  echo "[OK] FASTA index created: $REF_FAI"
fi

if [[ -f "$REF_DICT" ]]; then
  echo "[SKIP] Sequence dictionary already exists: $REF_DICT"
else
  echo "[RUN] Creating sequence dictionary (.dict) with GATK"
  $GATK /gatk/gatk CreateSequenceDictionary \
      -R "$REF_FA" \
      -O "$REF_DICT"
  echo "[OK] Sequence dictionary created: $REF_DICT"
fi

# -------------------------
# STEP 3: ClinVar
# -------------------------
log_step "[STEP 3] ClinVar (GRCh38 VCF + tabix index)"

mkdir_if_missing "${REF_ROOT}/ClinVar"
CLINVAR_VCF="${REF_ROOT}/ClinVar/clinvar.vcf.gz"
CLINVAR_TBI="${REF_ROOT}/ClinVar/clinvar.vcf.gz.tbi"

download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
  "$CLINVAR_VCF"

if [[ -f "$CLINVAR_TBI" ]]; then
  echo "[SKIP] ClinVar tabix index already exists: $CLINVAR_TBI"
else
  echo "[RUN] Indexing ClinVar with tabix"
  $TABIX -p vcf "$CLINVAR_VCF"
  echo "[OK] ClinVar indexed: $CLINVAR_TBI"
fi

# -------------------------
# STEP 4: CIViC
# -------------------------
log_step "[STEP 4] CIViC (nightly download)"

mkdir_if_missing "${REF_ROOT}/CIVIC"
CIVIC_FILE="${REF_ROOT}/CIVIC/nightly-VariantSummaries.tsv"

download_if_missing \
  "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv" \
  "$CIVIC_FILE"

# Capture version from download date since CIViC nightly has no embedded version
CIVIC_DOWNLOAD_DATE="$(date +%Y-%m-%d)"

# -------------------------
# STEP 5: REVEL
# -------------------------
log_step "[STEP 5] REVEL (zip download)"

mkdir_if_missing "${REF_ROOT}/REVEL"
REVEL_ZIP="${REF_ROOT}/REVEL/revel_all_chromosomes.zip"

download_if_missing \
  "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" \
  "$REVEL_ZIP"

# -------------------------
# STEP 6: SpliceAI
# -------------------------
log_step "[STEP 6] SpliceAI (SNV + INDEL VCFs)"

# SpliceAI is gated behind an Illumina BaseSpace login and cannot be
# downloaded automatically. Files must be manually downloaded and placed
# in SPLICEAI_STAGING_DIR before running this script. See comments at the end of this script
#
# Expected files:
#   spliceai_scores.raw.snv.hg38.vcf.gz
#   spliceai_scores.raw.indel.hg38.vcf.gz

mkdir_if_missing "${REF_ROOT}/SpliceAI"

SPLICEAI_SNV="${REF_ROOT}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="${REF_ROOT}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"

SPLICEAI_SNV_SRC="${SPLICEAI_STAGING_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL_SRC="${SPLICEAI_STAGING_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz"

SPLICEAI_MISSING=0

for src_file in "$SPLICEAI_SNV_SRC" "$SPLICEAI_INDEL_SRC"; do
  if [[ ! -f "$src_file" ]]; then
    echo "[WARN] SpliceAI staging file not found: $src_file"
    SPLICEAI_MISSING=1
  fi
done

if (( SPLICEAI_MISSING == 1 )); then
  echo ""
  echo "[WARN] =================================================================="
  echo "[WARN] SpliceAI files are missing from the staging directory:"
  echo "[WARN]   ${SPLICEAI_STAGING_DIR}"
  echo "[WARN]"
  echo "[WARN] SpliceAI requires a free Illumina BaseSpace account."
  echo "[WARN] Download both files manually from:"
  echo "[WARN]   SNV   : https://basespace.illumina.com/s/otSPW8hnhaZR"
  echo "[WARN]   INDEL : https://basespace.illumina.com/s/jI4aHqRKoGqE"
  echo "[WARN]"
  echo "[WARN] Then place them in: ${SPLICEAI_STAGING_DIR}"
  echo "[WARN] and re-run this script."
  echo "[WARN] =================================================================="
else
  declare -A SPLICEAI_FILES=(
    ["$SPLICEAI_SNV_SRC"]="$SPLICEAI_SNV"
    ["$SPLICEAI_INDEL_SRC"]="$SPLICEAI_INDEL"
  )

  for src in "${!SPLICEAI_FILES[@]}"; do
    dest="${SPLICEAI_FILES[$src]}"
    if [[ -f "$dest" ]]; then
      echo "[SKIP] Already exists: $dest"
    else
      echo "[CP] $src -> $dest"
      cp "$src" "$dest"
      echo "[OK] Copied: $dest"
    fi
  done
fi

# -------------------------
# STEP 7: gnomAD constraints (LOEUF)
# -------------------------
log_step "[STEP 7] gnomAD constraint metrics (LOEUF)"

mkdir_if_missing "${REF_ROOT}/gnomAD_constraints"

LOEUF_GRCH38="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_GRCH38_TBI="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"

echo "[INFO] Downloading gnomAD LOEUF from Google Drive."
echo "       NOTE: Google Drive may silently return an HTML error page instead"
echo "             of the file if the link has expired or the quota is exceeded."
echo "             Verify the downloaded file is not an HTML page if indexing fails."

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1pJMBwZQQLB2K4DelrWCEILzSadmLVybq" \
  "$LOEUF_GRCH38"

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1ZFg9d2VW_PC3VcsaSS84LaKeudw41PT8" \
  "$LOEUF_GRCH38_TBI"

# Sanity-check: the file should be binary (gzip), not HTML
if file "$LOEUF_GRCH38" | grep -q "HTML"; then
  echo "[ERROR] LOEUF download appears to be an HTML page, not a gzip file."
  echo "        The Google Drive link may have expired or hit a download quota."
  echo "        Remove ${LOEUF_GRCH38} and obtain the file manually, then re-run."
  rm -f "$LOEUF_GRCH38"
  exit 1
fi

# -------------------------
# STEP 8: VEP cache
# -------------------------
log_step "[STEP 8] VEP cache (homo_sapiens, GRCh38)"

VEP_CACHE_DIR="${REF_ROOT}/homo_sapiens"

if [[ -d "$VEP_CACHE_DIR" ]] && [[ -n "$(ls -A "$VEP_CACHE_DIR" 2>/dev/null)" ]]; then
  echo "[SKIP] VEP cache directory already present and non-empty: $VEP_CACHE_DIR"
else
  echo "[RUN] Installing VEP cache into: ${REF_ROOT}"
  echo "      (this may take a long time — cache is ~15 GB)"
  # Use the vep_install script available in the container's PATH rather
  # than a hardcoded miniconda path, which is version- and system-specific
  $VEP perl /opt/vep/src/ensembl-vep/INSTALL.pl \
    -a cf \
    -s homo_sapiens \
    -y GRCh38 \
    -c "${REF_ROOT}" \
    --NO_HTSLIB \
    --NO_TEST \
    --CACHE_VERSION 114
echo "[OK] VEP cache installed: $VEP_CACHE_DIR"
fi

# -------------------------
# STEP 9: VEP plugins
# -------------------------
log_step "[STEP 9] VEP Plugins"

PLUGINS_DIR="${REF_ROOT}/Plugins"

if [[ -d "${PLUGINS_DIR}/.git" ]]; then
  echo "[SKIP] VEP plugins repo already cloned: $PLUGINS_DIR"
else
  if [[ -d "$PLUGINS_DIR" ]]; then
    echo "[WARN] $PLUGINS_DIR exists but is not a git repo."
    echo "       Remove the directory manually if you want to re-clone, then re-run."
  else
    echo "[RUN] Cloning Ensembl VEP plugins repository..."
    git clone https://github.com/Ensembl/VEP_plugins.git "$PLUGINS_DIR"
    echo "[OK] VEP plugins cloned: $PLUGINS_DIR"
  fi
fi

# -------------------------
# STEP 10: Cancer Hotspots
# -------------------------
log_step "[STEP 10] Cancer Hotspots (hg38)"

mkdir_if_missing "${REF_ROOT}/CancerHotSpots"

HOTSPOT_VCF="${REF_ROOT}/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"

download_if_missing \
  "https://raw.githubusercontent.com/charlottekyng/cancer_hotspots/master/resources/hg38.hotspots_changv2_gao_nc.vcf.gz" \
  "$HOTSPOT_VCF"

# -------------------------
# STEP 11: Write config file
# -------------------------
log_step "[STEP 11] Writing pipeline configuration file"

CONFIG_FILE="pipeline_config.sh"

# Resource versions / metadata
INSTALL_DATE="$(date +%Y-%m-%d)"
REFERENCE_VERSION="UCSC_hg38"
CHAIN_VERSION="UCSC_hg19ToHg38"
CIVIC_VERSION="civic_nightly_${CIVIC_DOWNLOAD_DATE}"
REVEL_VERSION="v1.3"
SPLICEAI_VERSION="SpliceAIv1.3"
LOEUF_VERSION="custom_grch38"
CANCER_HOTSPOTS_VERSION="charlottekyng/cancer_hotspots_master"
VEP_CACHE_VERSION="GRCh38"
VEP_PLUGINS_VERSION="git_clone_$(date +%Y-%m-%d)"

CLINVAR_VERSION="$(zgrep -m1 '^##fileDate=' "$CLINVAR_VCF" 2>/dev/null | cut -d= -f2 || true)"
if [[ -z "$CLINVAR_VERSION" ]]; then
  CLINVAR_VERSION="unknown_check_${INSTALL_DATE}"
fi

cat > "$CONFIG_FILE" <<EOF
# --------------------------------------------------
# Clinical Genomics Pipeline configuration
# Generated automatically by install_references.sh
# --------------------------------------------------
#
# Installed by : $(whoami)
# Install date : ${INSTALL_DATE}
# Install host : $(hostname)
# --------------------------------------------------

INSTALL_DATE="${INSTALL_DATE}"
REF_ROOT="$(realpath "$REF_ROOT")"

REF_GENOME="\${REF_ROOT}/reference/hg38.fa"
REF_GENOME_VERSION="${REFERENCE_VERSION}"

CHAIN_FILE="\${REF_ROOT}/liftover/hg19ToHg38.over.chain"
CHAIN_FILE_VERSION="${CHAIN_VERSION}"

CLINVAR_VCF="\${REF_ROOT}/ClinVar/clinvar.vcf.gz"
CLINVAR_VERSION="${CLINVAR_VERSION}"

CIVIC_FILE="\${REF_ROOT}/CIVIC/nightly-VariantSummaries.tsv"
CIVIC_VERSION="${CIVIC_VERSION}"

REVEL_FILE="\${REF_ROOT}/REVEL/revel_all_chromosomes.zip"
REVEL_VERSION="${REVEL_VERSION}"

SPLICEAI_SNV="\${REF_ROOT}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="\${REF_ROOT}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
SPLICEAI_VERSION="${SPLICEAI_VERSION}"

LOEUF_FILE="\${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_FILE_INDEX="\${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"
LOEUF_VERSION="${LOEUF_VERSION}"

VEP_CACHE_DIR="\${REF_ROOT}/homo_sapiens"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION}"

VEP_PLUGINS_DIR="\${REF_ROOT}/Plugins"
VEP_PLUGINS_VERSION="${VEP_PLUGINS_VERSION}"

CANCER_HOTSPOTS_VCF="\${REF_ROOT}/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
CANCER_HOTSPOTS_VERSION="${CANCER_HOTSPOTS_VERSION}"
EOF

echo "[OK] Configuration file created: $CONFIG_FILE"

# -------------------------
# Summary
# -------------------------
echo ""
echo "================================================================================"
echo " DONE"
echo "================================================================================"
echo " Installed by   : $(whoami)"
echo " Install root   : $(realpath "$REF_ROOT")"
echo " Config file    : $(realpath "$CONFIG_FILE")"
echo ""
echo " If a step failed, re-run this script — completed steps will be skipped."
echo "================================================================================"


###############################################################################
# SpliceAI Download Script (HPC)
# This is what I did in case this helps
# Following commands do:
#   1. Installs BaseSpace CLI
#   2. Authenticates (interactive)
#   3. Lists project and dataset
#   4. Downloads dataset using tmux (safe for SSH disconnects)
#
# Usage:
#   chmod +x download_spliceai.sh
#   ./download_spliceai.sh
###############################################################################

# set -e

# ###############################################################################
# # CONFIGURATION
# ###############################################################################

# PROJECT_ID="66029966"
# DATASET_ID="ds.20a701bc58ab45b59de2576db79ac8d0"
# OUTPUT_DIR="${HOME}/spliceai_staging"
# TMUX_SESSION="spliceai_download"

# ###############################################################################
# # STEP 1: Install BaseSpace CLI (if not already installed)
# ###############################################################################

# echo "Installing BaseSpace CLI..."

# # Create local bin directory
# mkdir -p "${HOME}/bin"

# # Download CLI only if not already present
# if [ ! -f "${HOME}/bin/bs" ]; then
#     wget https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs \
#         -O "${HOME}/bin/bs"
#     chmod +x "${HOME}/bin/bs"
# fi

# # Add to PATH if needed
# if [[ ":$PATH:" != *":${HOME}/bin:"* ]]; then
#     export PATH="${HOME}/bin:$PATH"
#     echo 'export PATH="${HOME}/bin:$PATH"' >> "${HOME}/.bashrc"
# fi

# ###############################################################################
# # STEP 2: Authenticate (interactive step)
# ###############################################################################

# echo "Authenticating with BaseSpace..."
# echo "Follow the URL printed below in your browser, then press Enter."

# # This will prompt a URL → open it in your local browser
# bs auth || true

# ###############################################################################
# # STEP 3: Verify project and dataset (optional but useful)
# ###############################################################################

# echo "Checking project..."
# bs list projects | grep "${PROJECT_ID}" || echo "Project not found"

# echo "Checking datasets..."
# bs list datasets --project-id "${PROJECT_ID}" || true

# ###############################################################################
# # STEP 4: Download using tmux (prevents job interruption)
# ###############################################################################

# echo "Starting download in tmux session: ${TMUX_SESSION}"

# # Create output directory
# mkdir -p "${OUTPUT_DIR}"

# # If tmux session already exists → attach
# if tmux has-session -t "${TMUX_SESSION}" 2>/dev/null; then
#     echo "Session already exists. Attaching..."
#     tmux attach -t "${TMUX_SESSION}"
# else
#     # Start new tmux session in background
#     tmux new -d -s "${TMUX_SESSION}" "
#         echo 'Downloading SpliceAI dataset...';
#         bs dataset download \
#             --id ${DATASET_ID} \
#             -o ${OUTPUT_DIR};
#         echo 'Download complete.';
#         bash
#     "

#     echo "Download started in tmux."
#     echo "To monitor progress:"
#     echo "  tmux attach -t ${TMUX_SESSION}"
# fi
# Dont forget to tabix them! tabix -p vcf

# ###############################################################################
# # END
# ###############################################################################

# echo "Script finished."