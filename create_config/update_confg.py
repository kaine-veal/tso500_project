#!/usr/bin/env python3

import yaml
import gzip
import re

INPUT_YAML = "ConfigORIGINAL.yaml"
OUTPUT_YAML = "Config.yaml"
CLINVAR_VCF = "/mnt/data1/vep_test/ClinVar/clinvar.vcf.gz"
HOTSPOTS_VCF = "/mnt/data1/vep_test/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"


# ----------------------------
# Tier definitions
# ----------------------------

DROP_FIELDS = {
"CHROM",
"ONCOKB_highestSensitiveLevel",
"ONCOKB_highestResistanceLevel",
"ONCOKB_highestDiagnosticImplicationLevel",
"ONCOKB_highestPrognosticImplicationLevel",
"ONCOKB_geneExist",
"ONCOKB_GENE_EXIST",
"ONCOKB_VARIANT_EXIST",
"ONCOKB_variantExist",
"ONCOKB_alleleExist",
"ONCOKB_ONCOGENIC",
"ONCOKB_oncogenic",
"ONCOKB_EFFECT",
"ONCOKB_mutationEffect.knownEffect",
"ONCOKB_treatments",
"ONCOKB_diagnosticImplications",
"ONCOKB_prognosticImplications",
"ONCOKB_dataVersion",
"ONCOKB_hotspot",
"ONCOKB_geneSummary",
"ONCOKB_variantSummary",
"ONCOKB_tumorTypeSummary",
"ONCOKB_lastUpdate",
"ONCOKB_VUS.1",
"ONCOKB_query.referenceGenome",
"ONCOKB_query.hugoSymbol",
"ONCOKB_query.entrezGeneId",
"ONCOKB_query.alteration",
"ONCOKB_query.tumorType",
"CIViC_CSQ",
"ONCOKB_JSON",
"ClinVar_GENEINFO_",
}


TIER2_FIELDS = {
"NCBI_Build",
"Reference_Allele",
"Tumor_Seq_Allele1",
"Tumor_Seq_Allele2",
"CHROM",
"POS",
"ID",
"QUAL",
"FILTER",
"FORMAT",
"FORMAT_DATA",
"NM_Transcript",
"Allele",
"SYMBOL",
"Gene",
"Feature_type",
"Feature",
"BIOTYPE",
"CDS_position",
"Protein_position",
"Amino_acids",
"Codons",
"SYMBOL_SOURCE",
"HGNC_ID",
"CANONICAL",
"TSL",
"APPRIS",
"CCDS",
"ENSP",
"SWISSPROT",
"TREMBL",
"UNIPARC",
"UNIPROT_ISOFORM",
"SOURCE",
"GENE_PHENO",
"DOMAINS",
"miRNA",
"HGVS_OFFSET",
"gnomADe_AF",
"gnomADe_AFR_AF",
"gnomADe_AMR_AF",
"gnomADe_ASJ_AF",
"gnomADe_EAS_AF",
"gnomADe_FIN_AF",
"gnomADe_MID_AF",
"gnomADe_NFE_AF",
"gnomADe_REMAINING_AF",
"gnomADe_SAS_AF",
"gnomADg_AF",
"gnomADg_AFR_AF",
"gnomADg_AMI_AF",
"gnomADg_AMR_AF",
"gnomADg_ASJ_AF",
"gnomADg_EAS_AF",
"gnomADg_FIN_AF",
"gnomADg_MID_AF",
"gnomADg_NFE_AF",
"gnomADg_REMAINING_AF",
"gnomADg_SAS_AF",
"MAX_AF",
"MAX_AF_POPS",
"CLIN_SIG",
"SOMATIC",
"PHENO",
"PUBMED",
"MOTIF_NAME",
"MOTIF_POS",
"HIGH_INF_POS",
"MOTIF_SCORE_CHANGE",
"TRANSCRIPTION_FACTORS",
"SpliceAI_cutoff",
"SpliceAI_pred_DP_AG",
"SpliceAI_pred_DP_AL",
"SpliceAI_pred_DP_DG",
"SpliceAI_pred_DP_DL",
"SpliceAI_pred_SYMBOL",
"ClinVar",
"CIViC",
"CIViC_GN",
"CIViC_VT",
"CIViC_CSQ",
"CancerHotspots_HOTSPOT_GENE",
"CancerHotspots_HOTSPOT_HGVSp",
"CancerHotspots_HOTSPOT3D",
"CancerHotspots_HOTSPOT3D_GENE",
"CancerHotspots_HOTSPOT3D_HGVSp",
"CancerHotspots_HOTSPOTNC",
"CancerHotspots_HOTSPOTNC_GENE",
"CancerHotspots_HOTSPOTNC_HGVSc",
"BND_DEPTH",
"CIEND",
"CIGAR",
"CIPOS",
"DP",
"DUPSVLEN",
"END",
"FractionInformativeReads",
"HOMLEN",
"HOMSEQ",
"IMPRECISE",
"LEFT_SVINSSEQ",
"MATEID",
"MATE_BND_DEPTH",
"MQ",
"OriginalContig",
"OriginalStart",
"RIGHT_SVINSSEQ",
"ReverseComplementedAlleles",
"SVLEN",
"SVTYPE"
}


# --------------------------------
# Tier 3 placeholder
# --------------------------------
# Add fields here later when the team defines Tier 3
TIER3_FIELDS = {
# "Example_field",
# "Example_field2",
}


# ------------------------------------------------
# Extract ClinVar INFO descriptions from VCF header
# ------------------------------------------------

def load_clinvar_descriptions(vcf_path):

    descriptions = {}

    with gzip.open(vcf_path, "rt") as f:
        for line in f:

            if not line.startswith("##INFO="):
                continue

            # Extract INFO ID
            id_match = re.search(r'ID=([^,]+)', line)

            # Extract description
            desc_match = re.search(r'Description="([^"]+)"', line)

            if id_match and desc_match:
                info_id = id_match.group(1)
                description = desc_match.group(1)

                descriptions[info_id] = description

    print(f"Loaded {len(descriptions)} ClinVar INFO fields")

    return descriptions

# ------------------------------------------------
# Update YAML descriptions using ClinVar VCF header
# ------------------------------------------------

def update_clinvar_descriptions(fields, clinvar_info):

    updated = 0
    matched = 0

    for field_name, meta in fields.items():

        if not field_name.startswith("ClinVar_"):
            continue

        clinvar_key = field_name.replace("ClinVar_", "")

        if clinvar_key in clinvar_info:

            matched += 1
            new_description = clinvar_info[clinvar_key]

            if meta.get("description") != new_description:
                meta["description"] = new_description
                updated += 1

    print(f"ClinVar fields matched: {matched}")
    print(f"ClinVar descriptions updated: {updated}")


# ----------------------------
# Main script
# ----------------------------

def apply_tiers(input_yaml, output_yaml):

    with open(input_yaml) as f:
        data = yaml.safe_load(f)

    fields = data.get("fields", {})

    counts = {"tier1": 0, "tier2": 0, "tier3": 0, "drop": 0}

    for k, v in fields.items():
        if v is None:
            print(f"WARNING: field '{k}' has no metadata")

    for field, meta in fields.items():

        field_clean = field.strip()

        if field_clean in DROP_FIELDS:
            #print(field_clean)
            #print(DROP_FIELDS)
            meta["tier"] = "drop"
            counts["drop"] += 1
            continue

        if field_clean in TIER3_FIELDS:
            meta["tier"] = 3
            counts["tier3"] += 1

        elif field_clean in TIER2_FIELDS:
            meta["tier"] = 2
            counts["tier2"] += 1

        else:
            meta["tier"] = 1
            counts["tier1"] += 1


    # ----------------------------------
    # Update ClinVar descriptions
    # ----------------------------------

    clinvar_info = load_clinvar_descriptions(CLINVAR_VCF)

    update_clinvar_descriptions(fields, clinvar_info)


    # ----------------------------------
    # Write final YAML
    # ----------------------------------

    with open(output_yaml, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


    print("\nTier assignment complete\n")
    print("Summary:")
    print(f"Tier 1: {counts['tier1']}")
    print(f"Tier 2: {counts['tier2']}")
    print(f"Tier 3: {counts['tier3']}")
    print(f"Dropped: {counts['drop']}")


if __name__ == "__main__":
    apply_tiers(INPUT_YAML, OUTPUT_YAML)