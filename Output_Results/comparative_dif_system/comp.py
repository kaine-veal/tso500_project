#!/usr/bin/env python3
import argparse
import pandas as pd


KEY_COLS = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1"]


def read_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    # Normalize whitespace
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df


def variant_label(row: pd.Series) -> str:
    chrom = row["Chromosome"]
    pos = row["Start_Position"]
    ref = row["Reference_Allele"]
    alt = row["Tumor_Seq_Allele1"]
    return f"{chrom}:{pos}{ref}>{alt}"


def main():
    ap = argparse.ArgumentParser(
        description="Sort two TSVs by (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele1) and report diffs."
    )
    ap.add_argument("file1", help="First TSV")
    ap.add_argument("file2", help="Second TSV")
    ap.add_argument("--max", type=int, default=0, help="Max diff lines to print (0=all)")
    ap.add_argument("--ignore-cols", default="", help="Comma-separated columns to ignore in comparison")
    args = ap.parse_args()

    df1 = read_tsv(args.file1)
    df2 = read_tsv(args.file2)

    # 1) Check headers identical (names + order)
    if list(df1.columns) != list(df2.columns):
        raise SystemExit(
            "[ERROR] Headers differ (names and/or order). Cannot do strict comparison.\n"
            f"  file1 columns ({len(df1.columns)}): {list(df1.columns)}\n"
            f"  file2 columns ({len(df2.columns)}): {list(df2.columns)}"
        )

    # 2) Check key columns present
    missing = [c for c in KEY_COLS if c not in df1.columns]
    if missing:
        raise SystemExit(f"[ERROR] Missing required key columns: {missing}")

    # 3) (Optional) ensure same shape (you said they are same size)
    if df1.shape != df2.shape:
        raise SystemExit(f"[ERROR] Different shapes: file1={df1.shape}, file2={df2.shape}")

    ignore = [c.strip() for c in args.ignore_cols.split(",") if c.strip()]
    # Make sure we never ignore the key columns (it would break the label)
    ignore = [c for c in ignore if c not in KEY_COLS]

    compare_cols = [c for c in df1.columns if c not in ignore]

    # 4) Sort both by the key columns (stable)
    df1s = df1.sort_values(by=KEY_COLS, kind="mergesort").reset_index(drop=True)
    df2s = df2.sort_values(by=KEY_COLS, kind="mergesort").reset_index(drop=True)

    # 5) Compare cell-by-cell and print report
    diffs = 0
    printed = 0

    # Precompute variant labels from file1 after sort (assumes keys align after sorting)
    labels = df1s[KEY_COLS].apply(lambda r: f"{r['Chromosome']}:{r['Start_Position']}{r['Reference_Allele']}>{r['Tumor_Seq_Allele1']}", axis=1)

    for i in range(len(df1s)):
        # sanity check: key columns should match after sorting
        for kc in KEY_COLS:
            if df1s.iat[i, df1s.columns.get_loc(kc)] != df2s.iat[i, df2s.columns.get_loc(kc)]:
                raise SystemExit(
                    "[ERROR] After sorting, key columns do not align between files.\n"
                    f"Row {i+1} keys:\n"
                    f"  file1: {df1s.loc[i, KEY_COLS].to_dict()}\n"
                    f"  file2: {df2s.loc[i, KEY_COLS].to_dict()}\n"
                    "This usually means there are duplicate keys or different row sets."
                )

        label = labels.iat[i]

        # compare all columns (including key cols is pointless; they match by check above)
        for col in compare_cols:
            if col in KEY_COLS:
                continue
            v1 = df1s.at[i, col]
            v2 = df2s.at[i, col]
            if v1 != v2:
                diffs += 1
                if args.max == 0 or printed < args.max:
                    vv1 = v1 if v1 != "" else "<EMPTY>"
                    vv2 = v2 if v2 != "" else "<EMPTY>"
                    print(f"{label} ({col}) {vv1} --- {vv2}")
                    printed += 1

    # Total values compared (excluding key columns)
    compared_columns = [c for c in compare_cols if c not in KEY_COLS]
    total_compared = len(df1s) * len(compared_columns)

    print("\n--- SUMMARY ---")
    print(f"Rows: {len(df1s)}")
    print(f"Columns compared: {len(compared_columns)}")
    print(f"Total values compared: {total_compared}")
    print(f"Different values: {diffs} of {total_compared} compared")

    if args.max and diffs > args.max:
        print(f"Printed first {args.max} diffs (use --max 0 to print all)")

    raise SystemExit(0 if diffs == 0 else 1)


if __name__ == "__main__":
    main()