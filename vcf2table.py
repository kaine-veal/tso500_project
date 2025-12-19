#!/usr/bin/env python3

import csv
import re
import argparse


# ------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(
        description="Extract full VEP + OncoKB annotations for one selected transcript per variant"
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-t", "--transcripts", required=True, help="Transcript whitelist file (NM_*)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


# ------------------------------------------------------------
# Load NM transcript whitelist
# ------------------------------------------------------------
def load_nm_transcripts(filename):
    with open(filename) as f:
        return {line.strip().split(".")[0] for line in f if line.strip()}


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():

    args = get_args()

    # --------------------------------------------------------
    # Parse VCF header
    # --------------------------------------------------------
    CSQ_FIELDS = None
    VCF_COLUMNS = None

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: (.+)">', line)
                CSQ_FIELDS = m.group(1).split("|")

            elif line.startswith("#CHROM"):
                VCF_COLUMNS = line.strip().lstrip("#").split("\t")
                break

    if not CSQ_FIELDS:
        raise RuntimeError("CSQ header not found")

    nm_transcripts = load_nm_transcripts(args.transcripts)

    rows = []

    # --------------------------------------------------------
    # Parse variants
    # --------------------------------------------------------
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            values = line.rstrip("\n").split("\t")
            record = dict(zip(VCF_COLUMNS, values))

            info_raw = record.pop("INFO")

            # Parse INFO safely
            info = {}
            for part in info_raw.split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    info[k] = v
                else:
                    info[part] = True

            csq_raw = info.pop("CSQ", None)

            # Initialize row with all non-INFO VCF columns
            row = dict(record)

            # ------------------------------------------------
            # Handle CSQ
            # ------------------------------------------------
            selected_vep = {f: "" for f in CSQ_FIELDS}
            nm_selected = ""

            if csq_raw:
                for csq in csq_raw.split(","):
                    fields = csq.split("|")
                    if len(fields) < len(CSQ_FIELDS):
                        fields.extend([""] * (len(CSQ_FIELDS) - len(fields)))

                    vep = dict(zip(CSQ_FIELDS, fields))

                    nm = ""
                    for key in ("MANE_SELECT", "MANE", "MANE_PLUS_CLINICAL"):
                        if vep.get(key, "").startswith("NM_"):
                            nm = vep[key].split(".")[0]
                            break

                    if not nm and vep.get("HGVSc", "").startswith("NM_"):
                        nm = vep["HGVSc"].split(":")[0].split(".")[0]

                    if nm in nm_transcripts:
                        selected_vep = vep
                        nm_selected = nm
                        break

                if not nm_selected:
                    selected_vep = dict(zip(CSQ_FIELDS, fields))
                    nm_selected = nm

            row["NM_Transcript"] = nm_selected
            row.update(selected_vep)

            # ------------------------------------------------
            # Add remaining INFO fields (OncoKB, etc.)
            # ------------------------------------------------
            for k, v in info.items():
                row[k] = v

            rows.append(row)

    # --------------------------------------------------------
    # Build output columns
    # --------------------------------------------------------
    output_columns = (
        [c for c in VCF_COLUMNS if c != "INFO"] +
        ["NM_Transcript"] +
        CSQ_FIELDS +
        sorted({k for r in rows for k in r} -
               set(VCF_COLUMNS) -
               set(CSQ_FIELDS) -
               {"NM_Transcript"})
    )

    # --------------------------------------------------------
    # Write CSV
    # --------------------------------------------------------
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=output_columns)
        writer.writeheader()
        writer.writerows(rows)

    print(f"✓ Wrote {len(rows)} rows with {len(output_columns)} columns")


if __name__ == "__main__":
    main()
