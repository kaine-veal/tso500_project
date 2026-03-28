BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/FFG"

# MuTect2 VCFs
wget ${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz
wget ${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz
wget ${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz

# SomaticSniper VCFs
wget ${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz
wget ${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz
wget ${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz

# Strelka VCFs
wget ${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.strelka.vcf.gz
wget ${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.strelka.vcf.gz
wget ${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.strelka.vcf.gz


!/bin/bash
subsample_vcf.sh
For each .vcf.gz: decompress, delete the gz, then create a subsampled VCF
with 1-2 PASS variants per chromosome (~40 total) preserving the full header.

set -euo pipefail

VARIANTS_PER_CHR=2   # change to 1 if you want ~25 variants instead

for gz in *.vcf.gz; do

    base="${gz%.vcf.gz}"
    vcf="${base}.vcf"
    out="${base}.subsample.vcf"

    echo "────────────────────────────────────────"
    echo "Processing: $gz"

    # 1. Decompress
    echo "  Decompressing..."
    bgzip -d "$gz"          # produces ${base}.vcf
    # If bgzip is not available fall back to gunzip:
    # gunzip "$gz"

    # 2. Delete the gz
    # (bgzip -d already removes it; line below is a safety net if gunzip was used)
    [ -f "$gz" ] && rm -f "$gz"

    # 3. Extract header
    grep "^#" "$vcf" > "$out"

    # 4. Pick up to VARIANTS_PER_CHR PASS variants per chromosome
    #    - skip header lines (^#)
    #    - keep only lines where FILTER column (col 7) is exactly PASS
    #    - group by chromosome (col 1), take first N per group
    grep -v "^#" "$vcf" \
        | awk -v n="$VARIANTS_PER_CHR" '
            $7 == "PASS" {
                chr = $1
                count[chr]++
                if (count[chr] <= n) print
            }
        ' >> "$out"

    # 5. Count results
    total=$(grep -v "^#" "$out" | wc -l)
    chroms=$(grep -v "^#" "$out" | cut -f1 | sort -u | wc -l)
    echo "  Subsampled: $total variants across $chroms chromosomes → $out"

done

echo ""
echo "Done. Summary of subsampled files:"
echo "────────────────────────────────────────"
for out in *.subsample.vcf; do
    total=$(grep -v "^#" "$out" | wc -l)
    echo "  $out : $total variants"
done

for vcf in *.subsample.vcf; do     bgzip "$vcf" && bcftools index -t "${vcf}.gz";     echo "Done: ${vcf}.gz"; done