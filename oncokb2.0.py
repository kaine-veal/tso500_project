#!/usr/bin/env python3
"""
Combined OncoKB annotator for VEP-annotated VCF files (v2.4).

- SNVs / indels / MantaINS  → /annotate/mutations/byProteinChange
- CNVs (MantaDUP/DEL)       → /annotate/copyNumberAlterations

Requires:
    - ONCOKB_TOKEN environment variable (or hardcoded token)
    - VCF annotated by VEP with CSQ field
    - TSO500_transcripts_list.txt file in the same directory for transcript prioritization.
      (Matches NM IDs ignoring version numbers).
"""

import os
import re
import json
import sys
import argparse
import requests
from cyvcf2 import VCF, Writer

# --- FILE CONFIGURATION ---
PREFERRED_TRANSCRIPTS_FILE = "TSO500_transcripts_list.txt"

# VEP CSQ FIELD INDEXES (Based on the provided header)
MIN_FIELDS_FOR_MANE_STATUS = 42 
SYMBOL_INDEX = 3
FEATURE_INDEX = 6   # Ensembl Transcript ID (ENST...)

# --- REFSEQ/MANE INDICES ---
MANE_TAG_INDEX = 25              
MANE_SELECT_NM_INDEX = 26        
MANE_PLUS_CLINICAL_NM_INDEX = 27 
MANE_STATUS_INDEX = 41           

HGVSP_INDEX = 11    # Protein Change (p. notation)

# ---------------------------------------------------------
# 1. Map tumor names found in filenames → OncoTree codes
# ---------------------------------------------------------
CANCER_MAP = {
    "brain": "CNS",
    "breast": "BRCA",
    "cholangiocarcinoma": "CHOL",
    "colon": "COAD",
    "endometrial": "EMCA",
    "head-and-neck-squamous-cell-carcinoma-hnsc": "HNSC",
    "head-and-neck-snuc": "SNUC",
    "lung": "LUNG",
    "melanoma": "MEL",
    "ovarian": "OVT",
    "pancreatic": "PAAD",
    "prostate": "PRAD",
    "renal": "RCC",
    "sarcoma": "SARC",
    "thyroid": "THY",
    "thyroid-carcinoma": "THCA",
    "other": "OTHER",
}

ONCOKB_BASE = "https://www.oncokb.org/api/v1"
API_MUT_BY_PROTEIN = f"{ONCOKB_BASE}/annotate/mutations/byProteinChange"
API_CNA = f"{ONCOKB_BASE}/annotate/copyNumberAlterations"

REFERENCE_GENOME = "GRCh38"


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
ONCOKB_TOKEN = ""


def load_preferred_transcripts(file_path: str) -> set:
    """Load a set of preferred NM IDs from a file, stripping version numbers."""
    if not os.path.exists(file_path):
        sys.stderr.write(f"ERROR: Preferred transcripts file not found at {file_path}. Cannot prioritize transcripts.\n")
        return set()
    
    try:
        with open(file_path, 'r') as f:
            return {line.strip().split('.')[0] for line in f if line.strip() and not line.startswith('#')}
    except Exception as e:
        sys.stderr.write(f"ERROR reading preferred transcripts file: {e}\n")
        return set()

PREFERRED_TRANSCRIPTS = load_preferred_transcripts(PREFERRED_TRANSCRIPTS_FILE)


def clean_value(v):
    if v is None:
        return ""
    s = str(v)
    return s.replace("\n", " ").replace("\r", " ").replace(";", ",").strip()


def get_cancer_from_filename(filename: str) -> str:
    base = os.path.basename(filename)
    for part in re.split(r"[_.\-]+", base):
        key = part.lower()
        if key in CANCER_MAP:
            return CANCER_MAP[key]
    return "UNKNOWN"


def _extract_one_letter_protein(csq: str):
    fields = csq.split("|")
    try:
        protein = fields[HGVSP_INDEX]
    except IndexError:
        return None

    if not protein:
        return None

    if protein.startswith("p."):
        protein = protein[2:]

    protein = protein.replace("Ter", "*")

    aa3to1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
        "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U", "Pyl": "O",
    }

    m = re.match(r"([A-Za-z*]{3})(\d+)([A-Za-z*]{3})", protein)
    if m:
        ref3, pos, alt3 = m.groups()
        ref1 = aa3to1.get(ref3, ref3[0])
        alt1 = aa3to1.get(alt3, alt3[0])
        return f"{ref1}{pos}{alt1}"

    return protein


def find_preferred_csq_and_protein(variant, csqs: str, preferred_tx: set):
    """
    Finds the best CSQ annotation using a 3-tier prioritization logic.
    Returns: (gene, protein_change, selected_transcript, selection_reason)
    """
    if not csqs:
        return None, None, None, "No CSQ provided"

    annotations = csqs.split(',')
    selected_annotation = None
    selection_reason = "Fallback (VEP Default)"

    # --- TIER 1: PREFERRED TRANSCRIPTS (Match against NM IDs ignoring version) ---
    if preferred_tx:
        for ann in annotations:
            fields = ann.split('|')
            if len(fields) > MANE_SELECT_NM_INDEX:
                current_nm_id_select = fields[MANE_SELECT_NM_INDEX].split('.')[0].strip()
                if current_nm_id_select and current_nm_id_select in preferred_tx:
                    selected_annotation = ann
                    selection_reason = "Preferred List Match (NM/MANE Select)"
                    break

            if selected_annotation is None and len(fields) > MANE_PLUS_CLINICAL_NM_INDEX:
                current_nm_id_plus = fields[MANE_PLUS_CLINICAL_NM_INDEX].split('.')[0].strip()
                if current_nm_id_plus and current_nm_id_plus in preferred_tx:
                    selected_annotation = ann
                    selection_reason = "Preferred List Match (NM/MANE Plus Clinical)"
                    break
        if selected_annotation:
            pass 

    # --- TIER 2: MANE TRANSCRIPTS ---
    if selected_annotation is None and len(annotations) > 0 and len(annotations[0].split('|')) >= MIN_FIELDS_FOR_MANE_STATUS:
        for ann in annotations:
            fields = ann.split('|')
            if len(fields) > MANE_STATUS_INDEX and fields[MANE_STATUS_INDEX] == "MANE Select":
                selected_annotation = ann
                selection_reason = "MANE Select"
                break
        
        if selected_annotation is None:
            for ann in annotations:
                fields = ann.split('|')
                if len(fields) > MANE_STATUS_INDEX and fields[MANE_STATUS_INDEX].startswith("MANE Plus Clinical"):
                    selected_annotation = ann
                    selection_reason = "MANE Plus Clinical"
                    break
        
        if selected_annotation is None:
            selection_reason = "Fallback (VEP Default)"
            

    # --- TIER 3: FALLBACK ---
    if selected_annotation is None and annotations:
        selected_annotation = annotations[0]
        selection_reason = "Fallback (VEP Default)"

    if selected_annotation is None:
        selection_reason = "No valid CSQ found"
        return None, None, None, selection_reason

    fields = selected_annotation.split('|')
    gene = fields[SYMBOL_INDEX] if len(fields) > SYMBOL_INDEX else None
    selected_transcript = fields[FEATURE_INDEX] if len(fields) > FEATURE_INDEX else "N/A"
    protein_change = _extract_one_letter_protein(selected_annotation)
    
    sys.stdout.write(
        f"  [Variant {variant.CHROM}:{variant.POS} {variant.REF}>{variant.ALT}] "
        f"Selected: {selection_reason}. Gene: {gene}. Tx: {selected_transcript}. Protein: {protein_change}\n"
    )

    return gene, protein_change, selected_transcript, selection_reason


# ---------------------------------------------------------
# OncoKB queries
# ---------------------------------------------------------
def query_oncokb_mutation(gene: str, alteration: str, tumor_type: str,
                             token: str, cache: dict):
    if not gene or not alteration:
        return None

    key = ("mut", gene, alteration, tumor_type)
    if key in cache:
        return cache[key]

    params = {
        "referenceGenome": REFERENCE_GENOME,
        "hugoSymbol": gene,
        "alteration": alteration,
        "tumorType": tumor_type,
    }
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}

    resp = requests.get(API_MUT_BY_PROTEIN, params=params, headers=headers, timeout=15)
    if resp.status_code != 200:
        # Don't print 404s (not found), only actual errors
        if resp.status_code != 404:
            sys.stderr.write(
                f"OncoKB mutation API ERROR {resp.status_code} "
                f"for {gene} {alteration}, tumor={tumor_type}\n"
            )
        cache[key] = None
        return None

    try:
        data = resp.json()
    except json.JSONDecodeError:
        data = None

    cache[key] = data
    return data


def query_oncokb_cna(gene: str, cna_type: str, tumor_type: str,
                             token: str, cache: dict):
    if not gene or not cna_type:
        return None

    key = ("cna", gene, cna_type, tumor_type)
    if key in cache:
        return cache[key]

    params = {
        "referenceGenome": REFERENCE_GENOME,
        "hugoSymbol": gene,
        "copyNameAlterationType": cna_type,
        "tumorType": tumor_type,
    }
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}

    resp = requests.get(API_CNA, params=params, headers=headers, timeout=15)
    if resp.status_code != 200:
        if resp.status_code != 404:
            sys.stderr.write(
                f"OncoKB CNA API ERROR {resp.status_code} "
                f"for {gene} {cna_type}, tumor={tumor_type}\n"
            )
        cache[key] = None
        return None

    try:
        data = resp.json()
    except json.JSONDecodeError:
        data = None

    cache[key] = data
    return data


# ---------------------------------------------------------
# Flatten OncoKB JSON
# ---------------------------------------------------------
def summarize_oncokb(data: dict) -> dict:
    if not data:
        return {}

    out = {}
    out["oncogenic"] = clean_value(data.get("oncogenic"))
    out["hotspot"] = clean_value(data.get("hotspot"))
    out["vus"] = clean_value(data.get("vus"))

    out["gene_exist"] = clean_value(data.get("geneExist"))
    out["variant_exist"] = clean_value(data.get("variantExist"))

    out["dataVersion"] = clean_value(data.get("dataVersion"))
    out["lastUpdate"] = clean_value(data.get("lastUpdate"))

    me = data.get("mutationEffect") or {}
    out["effect_known"] = clean_value(me.get("knownEffect"))
    out["effect_desc"] = clean_value(me.get("description"))

    out["sens_level"] = clean_value(data.get("highestSensitiveLevel"))
    out["resist_level"] = clean_value(data.get("highestResistanceLevel"))
    out["diag_level"] = clean_value(data.get("highestDiagnosticImplicationLevel"))
    out["prog_level"] = clean_value(data.get("highestPrognosticImplicationLevel"))

    out["gene_summary"] = clean_value(data.get("geneSummary"))
    out["variant_summary"] = clean_value(data.get("variantSummary"))
    out["tumor_type_summary"] = clean_value(data.get("tumorTypeSummary"))
    out["diagnostic_summary"] = clean_value(data.get("diagnosticSummary"))
    out["prognostic_summary"] = clean_value(data.get("prognosticSummary"))

    treatments = data.get("treatments") or []
    short = []
    for t in treatments:
        level = t.get("level", "")
        drugs = t.get("drugs") or []
        names = "/".join(clean_value(d.get("drugName")) for d in drugs)
        short.append(f"{level}:{names}".strip(":"))

    out["treatments"] = "|".join([s for s in short if s])
    out["json"] = clean_value(json.dumps(data, separators=(",", ":")))

    return out


# ---------------------------------------------------------
# Classify variant type
# ---------------------------------------------------------
def classify_variant(variant) -> str | None:
    """
    Classifies variant based on INFO['SVTYPE'] or ID column (Column 3).
    Returns:
        "cnv": For MantaDUP, MantaDEL (Copy Number Alterations).
        "snv": For SNVs, INDELs, MantaINS, MantaBND (treated as mutations to check for protein change).
    """
    
    # 1. Try to get SVTYPE from INFO
    svtype = variant.INFO.get("SVTYPE")
    
    # 2. If valid SVTYPE not found or we want to be sure about Manta ID, check ID column
    variant_id = variant.ID if variant.ID else ""
    
    # Normalize SV type string
    sv_class = "snv" # Default is SNV (Mutation endpoint)

    if svtype:
        sv = str(svtype).upper()
        if sv in {"DUP", "DEL", "GAIN", "LOSS"}:
            sv_class = "cnv"
        elif sv in {"INS", "BND", "INV"}:
            sv_class = "snv" # Treat Structural insertions/fusions as mutations first
    
    # Fallback/Override based on Manta ID naming convention
    # This catches cases where SVTYPE might be missing or generic
    if "MantaDUP" in variant_id or "MantaDEL" in variant_id:
        sv_class = "cnv"
    elif "MantaINS" in variant_id or "MantaBND" in variant_id:
        sv_class = "snv"

    return sv_class


def map_svtype_to_cna(svtype: str, variant_id: str) -> str | None:
    """Map SVTYPE or ID → OncoKB copyNameAlterationType.
    OncoKB API: Does not understand "DUP" or "MantaDUP". 
    It has a strictly controlled vocabulary for Copy Number Alterations.
    It specifically expects the strings:
    AMPLIFICATION
    DELETION
    (It also accepts GAIN and LOSS, but AMPLIFICATION is the standard term for DUP).
    https://www.oncokb.org/swagger-ui/index.html#/Annotations/annotateCopyNumberAlterationsGetUsingGET_1
    """
    
    # Prioritize ID check for Manta consistency
    if variant_id:
        if "MantaDUP" in variant_id: return "AMPLIFICATION"
        if "MantaDEL" in variant_id: return "DELETION"

    if not svtype:
        return None
        
    sv = svtype.upper()
    if sv in {"DUP", "GAIN"}:
        return "AMPLIFICATION"
    if sv in {"DEL", "LOSS"}:
        return "DELETION"
        
    return None


# ---------------------------------------------------------
# MAIN FUNCTION — Annotate VCF
# ---------------------------------------------------------
def annotate_vcf(input_vcf: str, output_vcf: str, tumor_mode: str):
    token = ONCOKB_TOKEN

    if tumor_mode == 'filename':
        tumor_type = get_cancer_from_filename(input_vcf)
    elif tumor_mode == 'generic':
        tumor_type = "UNKNOWN" 
    else:
        sys.stderr.write(f"ERROR: Invalid tumor mode '{tumor_mode}'.\n")
        sys.exit(1)

    if tumor_type == "UNKNOWN" and tumor_mode == 'filename':
        sys.stderr.write(f"WARNING: Tumor type could not be determined from filename. Using '{tumor_type}'.\n")
    
    print(f"\n✔ Annotation Mode: {'Filename Inference' if tumor_mode == 'filename' else 'Generic Query'}")
    print(f"✔ Detected tumor type used for OncoKB: {tumor_type}")
    
    if PREFERRED_TRANSCRIPTS:
        print(f"✔ Using {len(PREFERRED_TRANSCRIPTS)} preferred NM transcripts (version agnostic).")
    else:
        print("⚠️ No preferred transcripts loaded.")

    vcf = VCF(input_vcf)

    fields_to_add = {
        "ONCOKB_ONCOGENIC": "OncoKB oncogenic classification",
        "ONCOKB_HOTSPOT": "OncoKB hotspot flag",
        "ONCOKB_EFFECT": "OncoKB mutationEffect.knownEffect",
        "ONCOKB_EFFECT_DESC": "OncoKB mutationEffect.description",
        "ONCOKB_SENS_LVL": "OncoKB highestSensitiveLevel",
        "ONCOKB_RESIST_LVL": "OncoKB highestResistanceLevel",
        "ONCOKB_DIAG_LVL": "OncoKB highestDiagnosticImplicationLevel",
        "ONCOKB_PROG_LVL": "OncoKB highestPrognosticImplicationLevel",
        "ONCOKB_TREATMENTS": "OncoKB treatment summary",
        "ONCOKB_GENE_SUMMARY": "OncoKB geneSummary",
        "ONCOKB_VARIANT_SUMMARY": "OncoKB variantSummary",
        "ONCOKB_TUMOR_TYPE_SUMMARY": "OncoKB tumorTypeSummary",
        "ONCOKB_DIAG_SUMMARY": "OncoKB diagnosticSummary",
        "ONCOKB_PROG_SUMMARY": "OncoKB prognosticSummary",
        "ONCOKB_VUS": "OncoKB VUS classification",
        "ONCOKB_GENE_EXIST": "OncoKB: Gene is curated in the database (boolean)",
        "ONCOKB_VARIANT_EXIST": "OncoKB: Exact variant/alteration is curated (boolean)",
        "ONCOKB_DATA_VERSION": "OncoKB data version",
        "ONCOKB_LAST_UPDATE": "OncoKB last update date",
        "ONCOKB_JSON": "Raw OncoKB JSON",
        "ONCOKB_QUERY_TYPE": "OncoKB query type: MUTATION or CNA",
    }

    for key, desc in fields_to_add.items():
        vcf.add_info_to_header({"ID": key, "Description": desc, "Type": "String", "Number": "1"})

    writer = Writer(output_vcf, vcf)

    total_variants = 0
    annotated_mut = 0
    annotated_cna = 0
    mut_cache = {}
    cna_cache = {}

    for variant in vcf:
        total_variants += 1
        
        # Determine Class based on SVTYPE and ID (Manta check)
        v_class = classify_variant(variant)

        csqs = variant.INFO.get("CSQ")
        gene = None
        protein_change = None

        if v_class == "snv":
            # Process as Mutation (SNV/Indel/INS/BND)
            gene, protein_change, selected_transcript, selection_reason = find_preferred_csq_and_protein(
                variant, str(csqs), PREFERRED_TRANSCRIPTS
            )

            if not gene or not protein_change:
                writer.write_record(variant)
                continue

            data = query_oncokb_mutation(gene, protein_change, tumor_type, token, mut_cache)
            summary = summarize_oncokb(data)

            if data:
                annotated_mut += 1
                variant.INFO["ONCOKB_QUERY_TYPE"] = "MUTATION"

        elif v_class == "cnv":
            # Process as CNA (DUP/DEL)
            if csqs:
                fields = str(csqs).split(',')[0].split('|')
                gene = fields[SYMBOL_INDEX] if len(fields) > SYMBOL_INDEX else None
            
            svtype = variant.INFO.get("SVTYPE")
            variant_id = variant.ID if variant.ID else ""
            cna_type = map_svtype_to_cna(str(svtype), variant_id)
            
            if not gene or not cna_type:
                writer.write_record(variant)
                continue

            data = query_oncokb_cna(gene, cna_type, tumor_type, token, cna_cache)
            summary = summarize_oncokb(data)

            if data:
                annotated_cna += 1
                variant.INFO["ONCOKB_QUERY_TYPE"] = "CNA"

        else:
            writer.write_record(variant)
            continue

        mapping = {
            "oncogenic": "ONCOKB_ONCOGENIC",
            "hotspot": "ONCOKB_HOTSPOT",
            "vus": "ONCOKB_VUS",
            "gene_exist": "ONCOKB_GENE_EXIST",     
            "variant_exist": "ONCOKB_VARIANT_EXIST", 
            "dataVersion": "ONCOKB_DATA_VERSION",
            "lastUpdate": "ONCOKB_LAST_UPDATE",
            "effect_known": "ONCOKB_EFFECT",
            "effect_desc": "ONCOKB_EFFECT_DESC",
            "sens_level": "ONCOKB_SENS_LVL",
            "resist_level": "ONCOKB_RESIST_LVL",
            "diag_level": "ONCOKB_DIAG_LVL",
            "prog_level": "ONCOKB_PROG_LVL",
            "gene_summary": "ONCOKB_GENE_SUMMARY",
            "variant_summary": "ONCOKB_VARIANT_SUMMARY",
            "tumor_type_summary": "ONCOKB_TUMOR_TYPE_SUMMARY",
            "diagnostic_summary": "ONCOKB_DIAG_SUMMARY",
            "prognostic_summary": "ONCOKB_PROG_SUMMARY",
            "treatments": "ONCOKB_TREATMENTS",
            "json": "ONCOKB_JSON",
        }

        for key, value in summary.items():
            if key in mapping and value:
                variant.INFO[mapping[key]] = value

        writer.write_record(variant)

    writer.close()

    print("\n📊 OncoKB Annotation Summary:")
    print(f"   Total variants processed:    {total_variants}")
    print(f"   Mutations annotated:         {annotated_mut}")
    print(f"   CNVs annotated:              {annotated_cna}")
    if total_variants:
        print(f"   Overall annotation rate:     {(annotated_mut + annotated_cna) / total_variants * 100:.2f}%")
    print(f"\n✔ OncoKB-annotated VCF written to: {output_vcf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate VCF files using OncoKB.")
    parser.add_argument("input_vcf", help="Input VCF file annotated by VEP.")
    parser.add_argument("output_vcf", help="Output VCF file.")
    parser.add_argument("--tumor_mode", choices=['filename', 'generic'], default='filename',
        help="Tumor type logic: 'filename' (infer from name) or 'generic' (UNKNOWN).")

    args = parser.parse_args()

    annotate_vcf(args.input_vcf, args.output_vcf, args.tumor_mode)