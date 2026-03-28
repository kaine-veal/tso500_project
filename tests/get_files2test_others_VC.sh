#!/bin/bash
# get_files2test_others_VC.sh
# Downloads VCF files, then subsamples ~1% of PASS variants
# from standard human chromosomes (chr1-22, chrX, chrY) only.

set -euo pipefail

BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/FFG"

# ---------------------------------------------------------------------------
# 1. Download
# ---------------------------------------------------------------------------

# MuTect2 VCFs
# curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"

# # SomaticSniper VCFs
# curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"

# # Strelka VCFs
# curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"

# ---------------------------------------------------------------------------
# 2. Decompress, filter and subsample each file
# ---------------------------------------------------------------------------

for gz in *.vcf.gz; do

    base="${gz%.vcf.gz}"
    vcf="${base}.vcf"
    out="${base}.subsample.vcf"

    echo "════════════════════════════════════════"
    echo "Processing: $gz"

    # Decompress
    echo "  Decompressing..."
    bgzip -d "$gz"
    [ -f "$gz" ] && rm -f "$gz"

    # Preserve full header
    grep "^#" "$vcf" > "$out"

    # Extract PASS variants on standard human chromosomes only
    # Accepted: chr1-chr22, chrX, chrY
    # Rejected: chrM, chrUn_*, *_random, *_alt, and any other non-standard contigs
    PASS_variants=$(grep -v "^#" "$vcf" \
        | awk '$7 == "PASS" && $1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/' )

    total=$(echo "$PASS_variants" | grep -c "." || true)

    if [[ "$total" -eq 0 ]]; then
        echo "  No PASS variants on standard chromosomes found — skipping."
        continue
    fi

    # Subsample ~1%: keep each variant with probability 0.01
    # Using awk rand() seeded per-file for reproducibility
    echo "$PASS_variants" \
        | awk -v seed="$RANDOM" 'BEGIN { srand(seed) } rand() < 0.01 { print }' \
        >> "$out"

    sampled=$(grep -v "^#" "$out" | wc -l)
    chroms=$(grep -v "^#" "$out" | cut -f1 | sort -u | tr '\n' ' ')
    echo "  Total PASS variants on standard chrs : $total"
    echo "  After 1% subsample                   : $sampled"
    echo "  Chromosomes present                  : $chroms"
    echo "  Output → $out"
    # Delete original full VCF to save space
    rm -f "$vcf"

done

echo ""
echo "════════════════════════════════════════"
echo "Done. Summary:"
for out in *.subsample.vcf; do
    total=$(grep -v "^#" "$out" | wc -l)
    echo "  $out : $total variants"
done
