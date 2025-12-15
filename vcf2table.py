#!/usr/bin/env python3

import csv
import re
import argparse


# ------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(description="Extract VEP+OncoKB annotations for selected NM transcripts (with fallback priority)")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-t", "--transcripts", required=True, help="Transcript list file (NM_...)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("--vep-fields", default="vep_fields.txt", help="File containing VEP annotation fields")
    parser.add_argument("--oncokb-fields", default="oncokb_fields.txt", help="File containing OncoKB annotation fields")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    return parser.parse_args()


# ------------------------------------------------------------
# Load list from file
# ------------------------------------------------------------
def load_list(filename):
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]


# ------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------
def main():

    args = get_args()

    print("\n=== extract_annotation.py ===")

    # -----------------------------
    # Load VEP + OncoKB field names
    # -----------------------------
    VEP_FIELDS = load_list(args.vep_fields)
    ONCOKB_FIELDS = load_list(args.oncokb_fields)

    print(f"Loaded {len(VEP_FIELDS)} VEP fields from: {args.vep_fields}")
    print(f"Loaded {len(ONCOKB_FIELDS)} OncoKB fields from: {args.oncokb_fields}")

    # -----------------------------
    # Load transcript whitelist
    # -----------------------------
    with open(args.transcripts) as f:
        nm_transcripts = {
            line.strip().split(".")[0]      # strip version
            for line in f if line.strip()
        }

    print(f"Loaded {len(nm_transcripts)} NM transcripts from: {args.transcripts}")

    rows = []

    # Debug counters
    total_variants = 0
    variants_with_csq = 0
    chosen_transcripts = 0
    fallback_used = 0

    # -----------------------------
    # Parse VCF
    # -----------------------------
    print(f"Reading VCF: {args.vcf}")

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            total_variants += 1
            fields = line.strip().split("\t")
            chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]

            # Extract VEP CSQ annotation
            m = re.search(r"CSQ=([^;]+)", info)
            if not m:
                # Variant has NO VEP annotation
                empty_vep = {field: "" for field in VEP_FIELDS}
                print(f"[NO CSQ] {chrom}:{pos} {ref}>{alt}")

                # Extract OncoKB fields if present
                onc = {}
                for key in ONCOKB_FIELDS:
                    m2 = re.search(rf"{key}=([^;]+)", info)
                    onc[key] = m2.group(1) if m2 else ""

                rows.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "NM_Transcript": "",
                    **empty_vep,
                    **onc
                })

                if args.debug:
                    print(f"[DEBUG] Variant without CSQ → included: {chrom}:{pos}")

                continue

            variants_with_csq += 1
            csq_entries = m.group(1).split(",")

            # Priority buckets
            best_entry = None
            mane_entries = []
            canonical_entries = []
            protein_coding_entries = []
            all_entries = []

            # -----------------------------
            # Evaluate all CSQ entries
            # -----------------------------
            for csq in csq_entries:
                csq_fields = csq.split("|")
                if len(csq_fields) < len(VEP_FIELDS):
                    continue

                vep = dict(zip(VEP_FIELDS, csq_fields))
                all_entries.append(vep)

                # Extract possible NM transcript
                nm_trans = ""

                for key in ("MANE_SELECT", "MANE", "MANE_PLUS_CLINICAL"):
                    val = vep.get(key, "")
                    if val.startswith("NM_"):
                        nm_trans = val.split(".")[0]
                        break

                if not nm_trans:
                    hgvsc = vep.get("HGVSc", "")
                    if hgvsc.startswith("NM_"):
                        nm_trans = hgvsc.split(":")[0].split(".")[0]

                # Case 1: NM transcript that matches whitelist → highest priority
                if nm_trans in nm_transcripts:
                    best_entry = (vep, nm_trans)
                    if args.debug:
                        print(f"[DEBUG] NM transcript match: {nm_trans}")
                    break

                # Fill fallback buckets
                if vep.get("MANE_SELECT", "").startswith("NM_"):
                    mane_entries.append(vep)

                if vep.get("CANONICAL") == "YES":
                    canonical_entries.append(vep)

                if vep.get("BIOTYPE") == "protein_coding":
                    protein_coding_entries.append(vep)

            # -----------------------------
            # Select best fallback entry
            # -----------------------------
            nm_trans = ""

            if not best_entry:
                fallback_used += 1

                if mane_entries:
                    vep = mane_entries[0]
                elif canonical_entries:
                    vep = canonical_entries[0]
                elif protein_coding_entries:
                    vep = protein_coding_entries[0]
                elif all_entries:
                    vep = all_entries[0]
                else:
                    continue  # shouldn't happen

                # Try assigning NM transcript if possible
                mane_val = vep.get("MANE_SELECT", "")
                if mane_val.startswith("NM_"):
                    nm_trans = mane_val.split(".")[0]

                best_entry = (vep, nm_trans)

                if args.debug:
                    print(f"[DEBUG] Fallback used: MANE/canonical/coding/first")

            # Final selected entry
            vep, nm_trans = best_entry
            chosen_transcripts += 1

            # -----------------------------
            # Extract OncoKB annotations
            # -----------------------------
            onc = {}
            for key in ONCOKB_FIELDS:
                m2 = re.search(rf"{key}=([^;]+)", info)
                onc[key] = m2.group(1) if m2 else ""

            # -----------------------------
            # Store final row
            # -----------------------------
            rows.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "NM_Transcript": nm_trans,
                **vep,
                **onc
            })

    # -----------------------------
    # Debug summary
    # -----------------------------
    print("\n=== SUMMARY ===")
    print(f"Total variants parsed:          {total_variants}")
    print(f"Variants with CSQ:              {variants_with_csq}")
    print(f"Selected transcript entries:    {chosen_transcripts}")
    print(f"Fallback selections used:       {fallback_used}")
    print(f"Rows to be written:             {len(rows)}")
    print("====================\n")

    # -----------------------------
    # Write output CSV
    # -----------------------------
    output_columns = (
        ["CHROM","POS","REF","ALT","NM_Transcript"] +
        VEP_FIELDS +
        ONCOKB_FIELDS
    )

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=output_columns)
        writer.writeheader()
        writer.writerows(rows)

    print(f"✓ Wrote {len(rows)} rows to {args.output}")
    print("✓ Done.\n")


# ------------------------------------------------------------
# ENTRY POINT
# ------------------------------------------------------------
if __name__ == "__main__":
    main()
