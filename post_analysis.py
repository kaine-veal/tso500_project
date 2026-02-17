import os
import pandas as pd

# To add a second header to inform the source of the annotation
HEADER_TO_SOURCE = {
    "NCBI_Build": "MAF Standard Metadata",
    "Hugo_Symbol": "MAF Standard Metadata",
    "Tumor_Sample_Barcode": "VCF",
    "HGVSp_Short": "VEP",
    "HGVSp": "VEP",
    "Chromosome": "VCF",
    "Start_Position": "VCF",
    "End_Position": "Computed",
    "Reference_Allele": "VCF",
    "Tumor_Seq_Allele1": "VCF",
    "Tumor_Seq_Allele2": "VCF",
    "CHROM": "VCF",
    "POS": "VCF",
    "ID": "VCF",
    "REF": "VCF",
    "ALT": "VCF",
    "QUAL": "VCF",
    "FILTER": "VCF",

    # VEP core annotation
    "NM_Transcript": "VEP",
    "Allele": "VEP",
    "Consequence": "VEP",
    "IMPACT": "VEP",
    "SYMBOL": "VEP",
    "Gene": "VEP",
    "Feature_type": "VEP",
    "Feature": "VEP",
    "BIOTYPE": "VEP",
    "EXON": "VEP",
    "INTRON": "VEP",
    "HGVSc": "VEP",
    "cDNA_position": "VEP",
    "CDS_position": "VEP",
    "Protein_position": "VEP",
    "Amino_acids": "VEP",
    "Codons": "VEP",
    "Existing_variation": "VEP",
    "DISTANCE": "VEP",
    "STRAND": "VEP",
    "FLAGS": "VEP",
    "VARIANT_CLASS": "VEP",
    "SYMBOL_SOURCE": "VEP",
    "HGNC_ID": "VEP",
    "CANONICAL": "VEP",
    "MANE": "VEP",
    "MANE_SELECT": "VEP",
    "MANE_PLUS_CLINICAL": "VEP",
    "TSL": "VEP",
    "APPRIS": "VEP",
    "CCDS": "VEP",
    "ENSP": "VEP",
    "SWISSPROT": "VEP",
    "TREMBL": "VEP",
    "UNIPARC": "VEP",
    "UNIPROT_ISOFORM": "VEP",
    "SOURCE": "VEP",
    "GENE_PHENO": "VEP",
    "SIFT": "VEP",
    "PolyPhen": "VEP",
    "DOMAINS": "VEP",
    "miRNA": "VEP",
    "HGVS_OFFSET": "VEP",
    "MOTIF_NAME": "VEP",
    "OriginalContig": "Manta",
    "OriginalStart": "Manta",

    # SpliceAI
    "SpliceAI_cutoff": "SpliceAI",
    "SpliceAI_pred_DS_AG": "SpliceAI",
    "SpliceAI_pred_DS_AL": "SpliceAI",
    "SpliceAI_pred_DS_DG": "SpliceAI",
    "SpliceAI_pred_DS_DL": "SpliceAI",
    "SpliceAI_pred_DP_AG": "SpliceAI",
    "SpliceAI_pred_DP_AL": "SpliceAI",
    "SpliceAI_pred_DP_DG": "SpliceAI",
    "SpliceAI_pred_DP_DL": "SpliceAI",
    "SpliceAI_pred_SYMBOL": "SpliceAI",


    "MOTIF_POS": "VEP",
    "HIGH_INF_POS": "VEP",
    "MOTIF_SCORE_CHANGE": "VEP",
    "TRANSCRIPTION_FACTORS": "VEP",

    "REVEL": "REVEL",
    "LOEUF": "gnomAD constraint metrics",

    # Population AFs
    "AF": "VEP",
    "AFR_AF": "gnomAD/1000G",
    "AMR_AF": "gnomAD/1000G",
    "EAS_AF": "gnomAD/1000G",
    "EUR_AF": "gnomAD/1000G",
    "SAS_AF": "gnomAD/1000G",

    # gnomAD exomes
    "gnomADe_AF": "gnomAD Exomes",
    "gnomADe_AFR_AF": "gnomAD Exomes",
    "gnomADe_AMR_AF": "gnomAD Exomes",
    "gnomADe_ASJ_AF": "gnomAD Exomes",
    "gnomADe_EAS_AF": "gnomAD Exomes",
    "gnomADe_FIN_AF": "gnomAD Exomes",
    "gnomADe_MID_AF": "gnomAD Exomes",
    "gnomADe_NFE_AF": "gnomAD Exomes",
    "gnomADe_REMAINING_AF": "gnomAD Exomes",
    "gnomADe_SAS_AF": "gnomAD Exomes",

    # gnomAD genomes
    "gnomADg_AF": "gnomAD Genomes",
    "gnomADg_AFR_AF": "gnomAD Genomes",
    "gnomADg_AMI_AF": "gnomAD Genomes",
    "gnomADg_AMR_AF": "gnomAD Genomes",
    "gnomADg_ASJ_AF": "gnomAD Genomes",
    "gnomADg_EAS_AF": "gnomAD Genomes",
    "gnomADg_FIN_AF": "gnomAD Genomes",
    "gnomADg_MID_AF": "gnomAD Genomes",
    "gnomADg_NFE_AF": "gnomAD Genomes",
    "gnomADg_REMAINING_AF": "gnomAD Genomes",
    "gnomADg_SAS_AF": "gnomAD Genomes",

    "MAX_AF": "VEP",
    "MAX_AF_POPS": "VEP",

    # ClinVar
    "CLIN_SIG": "ClinVar",
    "SOMATIC": "ClinVar",
    "PHENO": "ClinVar",
    "PUBMED": "ClinVar",
    "ClinVar": "ClinVar",
    "ClinVar_AF_ESP": "ClinVar",
    "ClinVar_AF_EXAC": "ClinVar",
    "ClinVar_AF_TGP": "ClinVar",
    "ClinVar_ALLELEID": "ClinVar",
    "ClinVar_CLNDN": "ClinVar",
    "ClinVar_CLNDNINCL": "ClinVar",
    "ClinVar_CLNDISDB": "ClinVar",
    "ClinVar_CLNDISDBINCL": "ClinVar",
    "ClinVar_CLNHGVS": "ClinVar",
    "ClinVar_CLNREVSTAT": "ClinVar",
    "ClinVar_CLNSIG": "ClinVar",
    "ClinVar_CLNSIGCONF": "ClinVar",
    "ClinVar_CLNSIGINCL": "ClinVar",
    "ClinVar_CLNSIGSCV": "ClinVar",
    "ClinVar_CLNVC": "ClinVar",
    "ClinVar_CLNVCSO": "ClinVar",
    "ClinVar_CLNVI": "ClinVar",
    "ClinVar_DBVARID": "ClinVar",
    "ClinVar_GENEINFO": "ClinVar",
    "ClinVar_MC": "ClinVar",
    "ClinVar_ONCDN": "ClinVar",
    "ClinVar_ONCDNINCL": "ClinVar",
    "ClinVar_ONCDISDB": "ClinVar",
    "ClinVar_ONCDISDBINCL": "ClinVar",
    "ClinVar_ONC": "ClinVar",
    "ClinVar_ONCINCL": "ClinVar",
    "ClinVar_ONCREVSTAT": "ClinVar",
    "ClinVar_ONCSCV": "ClinVar",
    "ClinVar_ONCCONF": "ClinVar",
    "ClinVar_ORIGIN": "ClinVar",
    "ClinVar_RS": "ClinVar",
    "ClinVar_SCIDN": "ClinVar",
    "ClinVar_SCIDNINCL": "ClinVar",
    "ClinVar_SCIDISDB": "ClinVar",
    "ClinVar_SCIDISDBINCL": "ClinVar",
    "ClinVar_SCIREVSTAT": "ClinVar",
    "ClinVar_SCI": "ClinVar",
    "ClinVar_SCIINCL": "ClinVar",
    "ClinVar_SCISCV": "ClinVar",

    # CIViC
    "CIViC": "CIViC",
    "CIViC_GN": "CIViC",
    "CIViC_VT": "CIViC",
    "CIViC_CSQ": "CIViC",

    # Cancer Hotspots
    "CancerHotspots": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT_GENE": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT_HGVSp": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT3D": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT3D_GENE": "Cancer Hotspots",
    "CancerHotspots_HOTSPOT3D_HGVSp": "Cancer Hotspots",
    "CancerHotspots_HOTSPOTNC": "Cancer Hotspots",
    "CancerHotspots_HOTSPOTNC_GENE": "Cancer Hotspots",
    "CancerHotspots_HOTSPOTNC_HGVSc": "Cancer Hotspots",
    "hotspot": "Cancer Hotspots",


    "ONCOKB_genefx": "OncoKB",
    "ONCOKB_variant": "OncoKB",
    "ONCOKB_allelef": "OncoKB",
    "ONCOKB_oncog": "OncoKB",
    "ONCOKB_highes": "OncoKB", 
    "ONCOKB_others": "OncoKB",  
    "ONCOKB_hotspot": "OncoKB",
    "ONCOKB_exon": "OncoKB",
    "ONCOKB_geneS": "OncoKB",  
    "ONCOKB_variantf": "OncoKB", 
    "ONCOKB_tumor1": "OncoKB",  
    "ONCOKB_progno": "OncoKB",  
    "ONCOKB_diagno": "OncoKB",  
    "ONCOKB_treatm": "OncoKB", 
    "ONCOKB_dataVe": "OncoKB",  
    "ONCOKB_lastUp": "OncoKB",
    "ONCOKB_highestSensitiveLevel": "OncoKB",
    "ONCOKB_highestResistanceLevel": "OncoKB",
    "ONCOKB_highestDiagnosticImplicationLevel": "OncoKB",
    "ONCOKB_highestPrognosticImplicationLevel": "OncoKB",
    "ONCOKB_highestFdaLevel": "OncoKB",
    "ONCOKB_exon": "OncoKB",
    "ONCOKB_prognosticSummary": "OncoKB",
    "ONCOKB_diagnosticSummary": "OncoKB",

    # Manta SV fields
    "BND_DEPTH": "Manta",
    "CIEND": "Manta",
    "CIGAR": "Manta",
    "CIPOS": "Manta",
    "DP": "VCF Standard",
    "DUPSVLEN": "Manta",
    "END": "Manta",
    "FractionInformativeReads": "Manta",
    "HOMLEN": "Manta",
    "HOMSEQ": "Manta",
    "IMPRECISE": "Manta",
    "LEFT_SVINSSEQ": "Manta",
    "MATEID": "Manta",
    "MATE_BND_DEPTH": "Manta",
    "MQ": "VCF Standard",
    "RIGHT_SVINSSEQ": "Manta",
    "ReverseComplementedAlleles": "Manta",
    "SVLEN": "Manta",
    "SVTYPE": "Manta",

    # OncoKB
    "ONCOKB_DATA_VERSION": "OncoKB",
    "ONCOKB_EFFECT": "OncoKB",
    "ONCOKB_GENE_EXIST": "OncoKB",
    "ONCOKB_GENE_SUMMARY": "OncoKB",
    "ONCOKB_HOTSPOT": "OncoKB",
    "ONCOKB_LAST_UPDATE": "OncoKB",
    "ONCOKB_ONCOGENIC": "OncoKB",
    "ONCOKB_QUERY_TYPE": "OncoKB",
    "ONCOKB_TUMOR_TYPE_SUMMARY": "OncoKB",
    "ONCOKB_VARIANT_EXIST": "OncoKB",
    "ONCOKB_VARIANT_SUMMARY": "OncoKB",
    "ONCOKB_treatments": "OncoKB",
    "ONCOKB_JSON": "OncoKB",

    # OncoKB query fields
    "ONCOKB_query.id": "OncoKB",
    "ONCOKB_query.referenceGenome": "OncoKB",
    "ONCOKB_query.hugoSymbol": "OncoKB",
    "ONCOKB_query.entrezGeneId": "OncoKB",
    "ONCOKB_query.alteration": "OncoKB",
    "ONCOKB_query.alterationType": "OncoKB",
    "ONCOKB_query.svType": "OncoKB",
    "ONCOKB_query.tumorType": "OncoKB",
    "ONCOKB_query.consequence": "OncoKB",
    "ONCOKB_query.proteinStart": "OncoKB",
    "ONCOKB_query.proteinEnd": "OncoKB",
    "ONCOKB_query.hgvs": "OncoKB",
    "ONCOKB_query.hgvsInfo": "OncoKB",
    "ONCOKB_query.canonicalTranscript": "OncoKB",

    # OncoKB mutation effect
    "ONCOKB_mutationEffect.knownEffect": "OncoKB",
    "ONCOKB_mutationEffect.description": "OncoKB",
    "ONCOKB_mutationEffect.citations.pmids": "OncoKB",
    "ONCOKB_mutationEffect.citations.abstracts": "OncoKB",

    # OncoKB summary fields
    "ANNOTATED": "Pipeline",
    "GENE_IN_ONCOKB": "OncoKB",
    "VARIANT_IN_ONCOKB": "OncoKB",
    "MUTATION_EFFECT": "OncoKB",
    "MUTATION_EFFECT_CITATIONS": "OncoKB",
    "ONCOGENIC": "OncoKB",
    "LEVEL_1": "OncoKB",
    "LEVEL_2": "OncoKB",
    "LEVEL_3A": "OncoKB",
    "LEVEL_3B": "OncoKB",
    "LEVEL_4": "OncoKB",
    "LEVEL_R1": "OncoKB",
    "LEVEL_R2": "OncoKB",
    "HIGHEST_LEVEL": "OncoKB",
    "HIGHEST_SENSITIVE_LEVEL": "OncoKB",
    "HIGHEST_RESISTANCE_LEVEL": "OncoKB",
    "TX_CITATIONS": "OncoKB",
    "LEVEL_Dx1": "OncoKB",
    "LEVEL_Dx2": "OncoKB",
    "LEVEL_Dx3": "OncoKB",
    "HIGHEST_DX_LEVEL": "OncoKB",
    "DX_CITATIONS": "OncoKB",
    "LEVEL_Px1": "OncoKB",
    "LEVEL_Px2": "OncoKB",
    "LEVEL_Px3": "OncoKB",
    "HIGHEST_PX_LEVEL": "OncoKB",
    "PX_CITATIONS": "OncoKB",


    # FORMAT + expanded fields
    "FORMAT": "VCF FORMAT",
    "FORMAT_DATA": "VCF FORMAT",
    "AD": "VCF FORMAT",
    "F1R2": "FORMAT",
    "F2R1": "FORMAT",
    "GT": "FORMAT",
    "MB": "FORMAT",
    "OBC": "FORMAT",
    "OBPa": "FORMAT",
    "OBParc": "FORMAT",
    "OBPsnp": "FORMAT",
    "PR": "FORMAT",
    "PS": "FORMAT",
    "SB": "FORMAT",
    "SQ": "FORMAT",
    "SR": "FORMAT",
    "VAF": "Pipeline",
}

# To filter out columns we don't need. On dev
cols_to_remove = ["CHROM", "FORMAT","FORMAT_DATA"]


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

# === Step 3: Calculate End_Position ===
def calculate_end_position(df):
    def compute_end(row):
        id_val = row.get("ID", "")
        ref = row.get("REF", "") or ""
        alt = row.get("ALT", "") or ""
        dupsvlen = row.get("DUPSVLEN", "")
        try:
            start = int(row["Start_Position"])
        except (ValueError, TypeError):
            return None

        if id_val == ".":
            if len(ref) == len(alt):
                return start
            return start + max(len(ref), len(alt)) - 1

        if id_val.startswith("MantaBND"):
            return None

        if id_val.startswith("MantaINS"):
            if alt == "<INS>":
                try:
                    return start + int(dupsvlen)
                except (ValueError, TypeError):
                    return None
            return start + len(alt)

        if id_val.startswith("MantaDUP"):
            try:
                return start + int(dupsvlen)
            except (ValueError, TypeError):
                return None

        return None

    df["End_Position"] = df.apply(compute_end, axis=1)
    return df

# === Step 3.1: Remove empty lists in OncoKB columns and the %3D issue ===
def clean_column_values(df):
    # Apply the simple string replacements across the whole DataFrame
    df = df.replace({"[]": "", "%3D": "="}, regex=False)

    # Transcript name included at start of cDNA change (HGVSc). Needs removing
    if "HGVSc" in df.columns:
        df["HGVSc"] = df["HGVSc"].astype(str).str.split(":").str[1]

    return df

# === Step 3.2: VAF from the AD columns ===
def add_vaf_column(df, ad_column="AD"):
    def compute_vaf(value):
        try:
            ref, alt = value.split(",")
            ref = float(ref)
            alt = float(alt)
            total = ref + alt
            return alt / total if total > 0 else 0
        except:
            return None  # or 0 if you prefer

    df["VAF"] = df[ad_column].apply(compute_vaf)
    return df



# === Step 4: Add source header row ===
def add_source_header_row(df, header_to_source):
    """
    To clarify where the data comes from

    This relies on the hard‑coded logic above, based on a few checks and quick searches

    It may still contain some errors
    """
    source_row = [header_to_source.get(col, "UNKNOWN") for col in df.columns]
    source_df = pd.DataFrame([source_row], columns=df.columns)
    df_with_source = pd.concat([source_df, df], ignore_index=True)
    return df_with_source

# === Step 4.1: Remove columns ===
def drop_columns(df, cols_to_remove): 
    return df.drop(columns=cols_to_remove, errors="ignore")

# === Step 5: Find empty columns in body ===
def find_empty_columns(df):
    """
    Look for columns whose body (rows 2 onward) is completely empty.
    Return list of (SOURCE, HEADER). This is just to show
    in the terminal completly empty column
    """
    empty_columns = []

    # Row 0 = SOURCE row
    source_row = df.iloc[0]

    # Body starts at row 2
    body = df.iloc[2:]

    for col in df.columns:
        col_body = body[col]

        # Treat NaN, empty strings, and whitespace as empty
        col_body_str = col_body.astype(str).str.strip()
        is_empty = col_body_str.replace({"": pd.NA, "nan": pd.NA}).isna().all()

        if is_empty:
            source = source_row[col]
            header = col  # the real header is the column name
            empty_columns.append((source, header))

    return empty_columns


# === Main merge function ===
def merge_maf_files(input_dir, output_file):
    maf_files = [f for f in os.listdir(input_dir) if f.lower().endswith(".maf")]
    if not maf_files:
        raise ValueError("No MAF files found in the directory.")

    dfs = []
    for maf in maf_files:
        path = os.path.join(input_dir, maf)
        print(f"Reading: {path}")
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
        df = standardize_format_column(df)
        dfs.append(df)

    merged = pd.concat(dfs, axis=0, ignore_index=True, sort=False)
    merged = calculate_end_position(merged)
    merged = normalize_format_fields(merged)
    merged = add_source_header_row(merged, HEADER_TO_SOURCE)
    merged = clean_column_values(merged)
    merged = add_vaf_column(merged)
    merged = drop_columns(merged,cols_to_remove )

    empties = find_empty_columns(merged)
    for source, header in empties:
        print(f"Empty column → SOURCE: {source} | HEADER: {header}")

    print("\n=== Final merged header (columns) ===")
    print("=====================================\n")
    print(f"Writing merged table to: {output_file}")
    merged.to_csv(output_file, sep="\t", index=False)


# === Entry point ===
if __name__ == "__main__":
    merge_maf_files(
        input_dir="./FINAL_Table/", #"./tmp/", # "./FINAL_Table/"
        output_file="merged_output.maf"
    )



