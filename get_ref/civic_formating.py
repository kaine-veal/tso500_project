"""
civic_tsv_to_vcf.py  —  Convert a CIViC nightly-VariantSummaries.tsv to a
                         bgzip-compressed, tabix-indexed VCF for use with VEP
                         --custom annotation.

Usage
-----
    python civic_tsv_to_vcf.py \\
        --input  nightly-VariantSummaries.tsv \\
        --output civic.vcf.gz \\
        [--assembly grch38]          # grch38 (default) or grch37

Requirements
------------
    Python >= 3.8  (standard library only; no third-party packages needed)
    bgzip and tabix must be on PATH (or passed via --bgzip / --tabix)

Output
------
    <output>          — bgzip-compressed VCF
    <output>.tbi      — tabix index

Column mapping (CIViC nightly TSV → VCF)
-----------------------------------------
GRCh38 mode (default):
    CHROM  ← chromosome (stripped of 'chr' prefix if present)
    POS    ← start38   (1-based; CIViC is already 1-based)
    REF    ← reference_bases  (. → N for variants without explicit alleles)
    ALT    ← variant_bases    (. → <CNV> or similar placeholder)

GRCh37 mode:
    POS    ← start      (GRCh37 coordinate)

All remaining columns are written as semicolon-separated key=value pairs
in the INFO field so VEP can pull them with --custom.

Notes on CIViC TSV format
--------------------------
The nightly release contains both GRCh37 and GRCh38 coordinates.
The exact column names have been stable since ~2020:

    variant_id, gene, entrez_id, variant, summary, variant_types,
    hgvs_expressions, last_review_date, civic_variant_evidence_score,
    allele_registry_id, clinvar_ids, variant_aliases,
    chromosome, start, stop, reference_bases, variant_bases,        ← GRCh37
    chromosome2, start2, stop2,                                      ← GRCh37 secondary (fusions)
    representative_transcript, ensembl_version, reference_build,
    chromosome_grch38, start_grch38, stop_grch38,                   ← GRCh38
    chromosome2_grch38, start2_grch38, stop2_grch38,                ← GRCh38 secondary (fusions)
    representative_transcript2, gene_civic_url, variant_civic_url

If the nightly file has changed column names, the script will print the
actual header and raise an informative error rather than silently producing
garbage output.
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
from datetime import date


# ---------------------------------------------------------------------------
# Accepted column name sets  (normalised to lowercase, underscores)
# We check both the new (grch38-suffixed) and old column variants.
# ---------------------------------------------------------------------------

CHROM_COLS_38  = ("chromosome_grch38", "chr_grch38", "chrom_grch38")
START_COLS_38  = ("start_grch38",)
STOP_COLS_38   = ("stop_grch38",)

CHROM_COLS_37  = ("chromosome", "chr", "chrom")
START_COLS_37  = ("start",)
STOP_COLS_37   = ("stop",)

REF_COLS       = ("reference_bases", "ref", "ref_bases")
ALT_COLS       = ("variant_bases",   "alt", "alt_bases", "var_bases")

# Fields to carry through into the INFO column (order matters for readability)
INFO_FIELDS = [
    "variant_id", "gene", "variant", "variant_types", "hgvs_expressions",
    "civic_variant_evidence_score", "allele_registry_id", "clinvar_ids",
    "reference_build", "representative_transcript",
    "gene_civic_url", "variant_civic_url",
]


def normalise_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "_", name.strip().lower())


def find_col(header_norm: list, candidates: tuple, label: str) -> int:
    for cand in candidates:
        if cand in header_norm:
            return header_norm.index(cand)
    raise ValueError(
        f"Cannot find required column '{label}'.\n"
        f"Tried: {candidates}\n"
        f"Available (normalised): {header_norm}"
    )


def sanitise_info_value(v: str) -> str:
    """Escape characters that would break VCF INFO syntax."""
    if not v or v in (".", "N/A", "NA", "None", ""):
        return "."
    # Remove semicolons (INFO field separator) and equals signs
    v = v.replace(";", "|").replace("=", ":").replace(" ", "_")
    # Truncate very long values (e.g. HGVS expressions) to keep VCF readable
    if len(v) > 500:
        v = v[:497] + "..."
    return v


def allele_ok(base: str) -> bool:
    """Return True if the base string is a valid VCF allele (not empty/dot)."""
    return bool(base) and base not in (".", "N/A", "NA", "None", "")


def convert(input_path: str, output_path: str, assembly: str,
            bgzip_bin: str, tabix_bin: str) -> None:

    assembly = assembly.lower().replace("-", "")
    use_38 = assembly == "grch38"

    # -----------------------------------------------------------------------
    # Pass 1: read entire TSV into memory (file is ~600 KB, easily fits RAM)
    # -----------------------------------------------------------------------
    opener = gzip.open if input_path.endswith(".gz") else open
    with opener(input_path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        raw_header = reader.fieldnames
        if raw_header is None:
            raise ValueError("Input file appears to be empty.")

        header_norm = [normalise_col(c) for c in raw_header]

        print(f"[INFO] Detected {len(raw_header)} columns in TSV header.")
        print(f"[INFO] Assembly mode: {assembly}")

        # Locate required columns
        # Fall back to GRCh37 columns if GRCh38-specific columns are absent
        # (CIViC nightly format dropped grch38-suffixed columns in 2025)
        if use_38:
            grch38_available = any(c in header_norm for c in CHROM_COLS_38)
            if not grch38_available:
                print("[WARN] GRCh38 columns not found in TSV — falling back to GRCh37 coordinates.")
                use_38 = False

        if use_38:
            chrom_idx = find_col(header_norm, CHROM_COLS_38, "chromosome_grch38")
            start_idx = find_col(header_norm, START_COLS_38, "start_grch38")
        else:
            chrom_idx = find_col(header_norm, CHROM_COLS_37, "chromosome")
            start_idx = find_col(header_norm, START_COLS_37, "start")

        try:
            ref_idx = find_col(header_norm, REF_COLS, "reference_bases")
        except ValueError:
            ref_idx = None

        try:
            alt_idx = find_col(header_norm, ALT_COLS, "variant_bases")
        except ValueError:
            alt_idx = None

        # Locate INFO columns (best-effort; missing ones are silently skipped)
        info_indices = {}
        for field in INFO_FIELDS:
            fn = normalise_col(field)
            if fn in header_norm:
                info_indices[field] = header_norm.index(fn)

        rows = list(reader)

    print(f"[INFO] Read {len(rows)} data rows from TSV.")

    # -----------------------------------------------------------------------
    # Pass 2: build VCF records
    # -----------------------------------------------------------------------
    vcf_records = []
    skipped = 0

    for row in rows:
        values = list(row.values())

        chrom = values[chrom_idx].strip() if values[chrom_idx] else ""
        start = values[start_idx].strip() if values[start_idx] else ""

        # Skip rows with missing coordinates
        if not chrom or chrom in (".", "N/A", "NA") or \
           not start or start in (".", "N/A", "NA") or not start.isdigit():
            skipped += 1
            continue

        # Strip 'chr' prefix — VEP GRCh38 cache uses bare chromosome numbers
        chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)

        pos = int(start)

        ref = values[ref_idx].strip().upper() if ref_idx is not None else ""
        alt = values[alt_idx].strip().upper() if alt_idx is not None else ""

        # Handle missing alleles (CNVs, expression events, fusions, etc.)
        if not allele_ok(ref):
            ref = "N"
        if not allele_ok(alt):
            alt = f"<{values[info_indices.get('variant_types', 0)] or 'UNK'}>"
            alt = re.sub(r"[^A-Za-z0-9_<>]", "_", alt)[:20]  # safe placeholder
            if not alt.startswith("<"):
                alt = "<UNK>"

        # Build INFO string
        info_parts = []
        for field, idx in info_indices.items():
            val = sanitise_info_value(values[idx] if idx < len(values) else "")
            info_parts.append(f"CIVIC_{field.upper()}={val}")

        info_str = ";".join(info_parts) if info_parts else "."

        vcf_records.append((chrom, pos, ref, alt, info_str))

    print(f"[INFO] {len(vcf_records)} records will be written "
          f"({skipped} skipped — missing {assembly} coordinates).")

    # -----------------------------------------------------------------------
    # Sort by chromosome then position
    # -----------------------------------------------------------------------
    def chrom_sort_key(rec):
        chrom = rec[0]
        # Numeric chromosomes sort as integers; X=23, Y=24, M/MT=25
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        if chrom in mapping:
            return (mapping[chrom], rec[1])
        try:
            return (int(chrom), rec[1])
        except ValueError:
            return (99, rec[1])

    vcf_records.sort(key=chrom_sort_key)

    # -----------------------------------------------------------------------
    # Write plain VCF to a temp file, then bgzip + tabix
    # -----------------------------------------------------------------------
    tmp_vcf = output_path.replace(".gz", "") + ".tmp.vcf"
    today = date.today().isoformat()

    with open(tmp_vcf, "w") as out:
        # VCF header
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={today}\n")
        out.write(f"##source=CIViC_nightly_VariantSummaries\n")
        out.write(f"##reference={assembly}\n")

        # INFO meta-lines
        for field in INFO_FIELDS:
            key = f"CIVIC_{field.upper()}"
            out.write(
                f'##INFO=<ID={key},Number=1,Type=String,'
                f'Description="CIViC field: {field}">\n'
            )

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for chrom, pos, ref, alt, info in vcf_records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    print(f"[INFO] Plain VCF written to: {tmp_vcf}")

    # bgzip
    final_vcf = output_path if output_path.endswith(".gz") else output_path + ".gz"
    print(f"[RUN] bgzip -> {final_vcf}")
    subprocess.run([bgzip_bin, "-f", "-c", tmp_vcf],
                   stdout=open(final_vcf, "wb"), check=True)
    os.remove(tmp_vcf)

    # tabix
    print(f"[RUN] tabix -> {final_vcf}.tbi")
    subprocess.run([tabix_bin, "-p", "vcf", final_vcf], check=True)

    print(f"[DONE] Output: {final_vcf}")
    print(f"[DONE] Index:  {final_vcf}.tbi")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CIViC nightly-VariantSummaries.tsv to a VCF for VEP --custom"
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="Path to nightly-VariantSummaries.tsv (plain or .gz)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output path (e.g. civic.vcf.gz)")
    parser.add_argument("-a", "--assembly", default="grch38",
                        choices=["grch38", "grch37"],
                        help="Reference assembly to use for coordinates (default: grch38)")
    parser.add_argument("--bgzip", default="bgzip",
                        help="Path to bgzip binary (default: bgzip from PATH)")
    parser.add_argument("--tabix", default="tabix",
                        help="Path to tabix binary (default: tabix from PATH)")
    args = parser.parse_args()

    convert(
        input_path=args.input,
        output_path=args.output,
        assembly=args.assembly,
        bgzip_bin=args.bgzip,
        tabix_bin=args.tabix,
    )


if __name__ == "__main__":
    main()
