#!/usr/bin/env python3

import os
import yaml
import pandas as pd


COLS_TO_REMOVE = []


# === Load YAML config ===
def load_field_config(yaml_path):

    with open(yaml_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if not config or "fields" not in config:
        raise ValueError(f"YAML file '{yaml_path}' does not contain a top-level 'fields' section.")

    return config["fields"]


# === Step 1: Standardize FORMAT_DATA column per file ===
def standardize_format_column(df):

    if "FORMAT" not in df.columns:
        return df

    format_idx = df.columns.get_loc("FORMAT")

    if format_idx + 1 < len(df.columns):
        sample_col = df.columns[format_idx + 1]
        df = df.rename(columns={sample_col: "FORMAT_DATA"})

    return df


# === Step 2: Expand FORMAT fields ===
def normalize_format_fields(df):

    if "FORMAT" not in df.columns or "FORMAT_DATA" not in df.columns:
        return df

    def parse_format(row):

        format_str = str(row.get("FORMAT", "") or "")
        data_str = str(row.get("FORMAT_DATA", "") or "")

        keys = format_str.split(":")
        values = data_str.split(":")

        if len(values) < len(keys):
            values += [None] * (len(keys) - len(values))

        return dict(zip(keys, values))

    parsed_rows = df.apply(parse_format, axis=1)

    all_keys = set()

    for d in parsed_rows:
        all_keys.update(d.keys())

    for key in sorted(all_keys):

        if key in df.columns:
            continue

        df[key] = parsed_rows.apply(lambda d: d.get(key, None))

    return df


# === Step 2.5: Expand structured annotation columns (e.g CIViC_CSQ) ===
def expand_structured_columns(df, field_config):

    for field, meta in field_config.items():

        if field not in df.columns:
            continue

        description = meta.get("description", "")

        if "|" not in description:
            continue

        # Avoid splitting natural language descriptions
        tokens = description.split("|")

        if any(len(t.strip()) > 40 for t in tokens):
            continue

        subfields = [x.strip() for x in description.split("|")]
        n_fields = len(subfields)

        print("\n====================================================")
        print(f"Expanding structured column: {field}")
        print(f"Detected {n_fields} subfields")

        print("\nHeader structure:")
        for i, sf in enumerate(subfields):
            print(f"  {i:02d} -> {sf}")

        def split_first_record(value):

            if pd.isna(value):
                return [None] * n_fields

            value = str(value)

            # CASE 1 — VEP style (records separated by comma)
            if "," in value:
                first_record = value.split(",")[0]
                parts = first_record.split("&")

            # CASE 2 — CIViC style (records concatenated with &)
            else:
                parts_all = value.split("&")

                if len(parts_all) % n_fields != 0:
                    print("\nWARNING: Field count mismatch!")
                    print(f"Total values: {len(parts_all)} | Expected multiple of {n_fields}")

                records = [
                    parts_all[i:i+n_fields]
                    for i in range(0, len(parts_all), n_fields)
                ]

                parts = records[0]

            if len(parts) < n_fields:
                parts += [None] * (n_fields - len(parts))

            parts = parts[:n_fields]

            # Debug print
            print("\nValue check:")
            print("Raw value:", value[:200], "...")

            for i, (header, val) in enumerate(zip(subfields, parts)):
                print(f"  {i:02d} | {header:35} -> {val}")

            return parts

        expanded = df[field].apply(split_first_record).apply(pd.Series)

        expanded.columns = [f"{field}_{c}" for c in subfields]

        print("\nExpanded column names:")
        for c in expanded.columns:
            print(" ", c)

        df = df.drop(columns=[field])
        df = pd.concat([df, expanded], axis=1)

    return df

# === Step 3: Calculate End_Position ===
def calculate_end_position(df):

    def compute_end(row):

        id_val = row.get("ID", "")
        ref = row.get("REF", "") or ""
        alt = row.get("ALT", "") or ""
        dupsvlen = row.get("DUPSVLEN", "")

        try:
            start = int(row["Start_Position"])
        except (ValueError, TypeError, KeyError):
            return None

        if id_val == ".":
            if len(ref) == len(alt):
                return start
            return start + max(len(ref), len(alt)) - 1

        if str(id_val).startswith("MantaBND"):
            return None

        if str(id_val).startswith("MantaINS"):
            if alt == "<INS>":
                try:
                    return start + int(dupsvlen)
                except (ValueError, TypeError):
                    return None
            return start + len(alt)

        if str(id_val).startswith("MantaDUP"):
            try:
                return start + int(dupsvlen)
            except (ValueError, TypeError):
                return None

        return None

    df["End_Position"] = df.apply(compute_end, axis=1)

    return df


# === Step 4: Clean values ===
def clean_column_values(df):

    df = df.replace({"[]": "", "%3D": "="}, regex=False)

    if "HGVSc" in df.columns:
        df["HGVSc"] = (
            df["HGVSc"]
            .astype(str)
            .apply(lambda x: x.split(":", 1)[1] if ":" in x else x)
        )

    return df


# === Step 5: Compute VAF ===
def add_vaf_column(df, ad_column="AD"):

    if ad_column not in df.columns:
        return df

    def compute_vaf(value):

        try:
            ref, alt = str(value).split(",")
            ref = float(ref)
            alt = float(alt)
            total = ref + alt
            return alt / total if total > 0 else 0
        except Exception:
            return None

    df["VAF"] = df[ad_column].apply(compute_vaf)

    return df


# === Step 6: Drop columns ===
def drop_columns(df, cols_to_remove):

    return df.drop(columns=cols_to_remove, errors="ignore")


# === Step 7: Report undocumented columns ===
def report_undocumented_columns(df, field_config):

    documented = set(field_config.keys())
    observed = set(df.columns)

    undocumented = sorted(observed - documented)

    if undocumented:

        print("\n=== WARNING: Undocumented columns found ===")

        for col in undocumented:
            print(f"Undocumented column: {col}")


# === Step 8: Build two-row header ===
def build_multiindex_header(df, field_config):

    metadata_row = []

    for col in df.columns:

        meta = field_config.get(col, {})

        description = meta.get("description", "UNKNOWN")
        source = meta.get("source", "UNKNOWN")
        version = meta.get("version", "UNKNOWN")

        metadata = f"{description} | {source} | {version}"

        metadata_row.append(metadata)

    df.columns = pd.MultiIndex.from_arrays(
        [df.columns, metadata_row],
        names=["Field", "Metadata"]
    )

    return df


# === Step 9: Split dataframe by tier ===
def split_by_tier(df, field_config):

    tier1_cols = []
    tier2_cols = []
    drop_cols = []

    for col in df.columns:

        meta = field_config.get(col, {})
        tier = meta.get("tier")

        if tier == "drop":
            drop_cols.append(col)

        elif tier == 1:
            tier1_cols.append(col)
            tier2_cols.append(col)

        elif tier == 2:
            tier2_cols.append(col)

        else:
            tier1_cols.append(col)
            tier2_cols.append(col)

    print(f"\nTier1 columns: {len(tier1_cols)}")
    print(f"Tier2 columns: {len(tier2_cols)}")
    print(f"Dropped columns: {len(drop_cols)}")

    df_tier1 = df[tier1_cols]
    df_tier2 = df[tier2_cols]

    return df_tier1, df_tier2


# === Main merge function ===
def merge_maf_files(input_dir, output_file, yaml_config):

    field_config = load_field_config(yaml_config)

    maf_files = [f for f in os.listdir(input_dir) if f.lower().endswith(".maf")]

    if not maf_files:
        raise ValueError("No MAF files found in the directory.")

    dfs = []

    for maf in maf_files:

        path = os.path.join(input_dir, maf)

        print(f"Reading: {path}")

        df = pd.read_csv(path, sep="\t", dtype=str)

        df = standardize_format_column(df)

        dfs.append(df)

    merged = pd.concat(dfs, axis=0, ignore_index=True, sort=False)

    merged = calculate_end_position(merged)
    merged = normalize_format_fields(merged)
    merged = expand_structured_columns(merged, field_config)
    merged = clean_column_values(merged)
    merged = add_vaf_column(merged)
    merged = drop_columns(merged, COLS_TO_REMOVE)

    report_undocumented_columns(merged, field_config)

    tier1_df, tier2_df = split_by_tier(merged, field_config)

    merged = build_multiindex_header(merged, field_config)
    tier1_df = build_multiindex_header(tier1_df, field_config)
    tier2_df = build_multiindex_header(tier2_df, field_config)

    master_file = output_file
    tier1_file = output_file.replace(".maf", "_tier1.maf")
    tier2_file = output_file.replace(".maf", "_tier2.maf")

    print(f"\nWriting master table: {master_file}")
    merged.to_csv(master_file, sep="\t", index=False)
    

    print(f"Writing Tier1 table: {tier1_file}")
    tier1_df.to_csv(tier1_file, sep="\t", index=False)

    print(f"Writing Tier2 table: {tier2_file}")
    tier2_df.to_csv(tier2_file, sep="\t", index=False)


# === Entry point ===
if __name__ == "__main__":

    merge_maf_files(
        input_dir="./FINAL_Table/",
        output_file="merged_output.maf",
        yaml_config="/mnt/data1/vep_sample_test/manuel/DoItAll/create_config/Config.yaml"
    )