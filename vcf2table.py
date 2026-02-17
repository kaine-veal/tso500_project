#!/usr/bin/env python3

import csv
import re
import argparse
import pandas as pd


# ------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(
        description="Extract full VEP + OncoKB annotations for one selected transcript per variant"
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-t", "--transcripts", required=True, help="Transcript whitelist file (NM_*)")
    parser.add_argument("-o", "--output", required=True, help="Final output file")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


# ------------------------------------------------------------
# Load NM transcript whitelist
# ------------------------------------------------------------
def load_nm_transcripts(filename):
    with open(filename) as f:
        return {line.strip().split(".")[0] for line in f if line.strip()}


# ------------------------------------------------------------
# MAIN: parse VCF → return DataFrame
# ------------------------------------------------------------
def main(args):

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

            # Add remaining INFO fields
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

    # First: build the DataFrame
    df = pd.DataFrame(rows, columns=output_columns)

    # Nothing else here — keep FORMAT and sample column as-is
    return df




# ------------------------------------------------------------
# Additional requirements
# ------------------------------------------------------------
def calculate_end(row):
    pos = int(row["POS"])
    ref = str(row["REF"])
    alt = str(row["ALT"])
    variant_class = str(row.get("VARIANT_CLASS", "")).lower()

    if variant_class == "snv":
        return pos

    if variant_class == "chromosome_breakpoint":
        return ""

    if variant_class == "insertion":
        svlen = row.get("SVLEN")
        if pd.notna(svlen):
            try:
                svlen = int(str(svlen).split(",")[0])
                return pos + svlen - 1
            except ValueError:
                pass
        return pos + len(alt) - 1

    if variant_class == "deletion":
        return pos + len(ref) - 1

    return pos


aa_dict = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

def convert_hgvsp_short(hgvsp):
    if not isinstance(hgvsp, str) or "p." not in hgvsp:
        return ""
    match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z\*]+)", hgvsp)
    if match:
        ref, pos, alt = match.groups()
        return f"p.{aa_dict.get(ref, ref)}{pos}{aa_dict.get(alt, alt)}"
    return hgvsp.split(":")[-1]


def final_requiretments(df, output_csv):

    if "END" not in df.columns:
        df["END"] = df.apply(calculate_end, axis=1)

    fixed_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    sample_candidates = [c for c in df.columns if c not in fixed_cols]
    sample_col = sample_candidates[0] if sample_candidates else None

    df["NCBI_Build"] = "GRCh38"
    df["Tumor_Sample_Barcode"] = sample_col if sample_col else ""

    if "HGVSp" in df.columns:
        df["HGVSp_Short"] = df["HGVSp"].apply(convert_hgvsp_short)
    else:
        df["HGVSp_Short"] = ""

    df["Hugo_Symbol"] = df["SYMBOL"] if "SYMBOL" in df.columns else ""

    df["Chromosome"] = df["CHROM"]
    df["Start_Position"] = df["POS"]
    df["End_Position"] = df["END"]
    df["Reference_Allele"] = df["REF"]
    df["Tumor_Seq_Allele1"] = df["ALT"]
    df["Tumor_Seq_Allele2"] = df["ALT"]

    first_cols = [
        "NCBI_Build", "Hugo_Symbol", "Tumor_Sample_Barcode",
        "HGVSp_Short", "HGVSp", "Chromosome", "Start_Position",
        "End_Position", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
    ]

    first_cols = [c for c in first_cols if c in df.columns]
    remaining = [c for c in df.columns if c not in first_cols]

    df = df[first_cols + remaining]

    # --------------------------------------------------------
    # Expand ONCOKB_JSON into separate columns and move to end
    # --------------------------------------------------------
    if "ONCOKB_JSON" in df.columns:
        import json
        from pandas import json_normalize

        # Parse JSON safely
        parsed = df["ONCOKB_JSON"].apply(
            lambda x: json.loads(x) if isinstance(x, str) and x.strip().startswith("{") else {}
        )

        # Deep flatten JSON (max_level=None allows full expansion)
        expanded = json_normalize(parsed, sep=".", max_level=None)

        # Prefix to avoid collisions
        expanded.columns = [f"ONCOKB_{c}" for c in expanded.columns]

        # Rebuild dataframe: original columns → expanded JSON → raw JSON at end
        df_no_json = df.drop(columns=["ONCOKB_JSON"])
        df = pd.concat([df_no_json, expanded, df["ONCOKB_JSON"]], axis=1)
    print(output_csv)
    df.to_csv(output_csv, index=False, sep="\t")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    args = get_args()
    df = main(args)
    final_requiretments(df, args.output)
