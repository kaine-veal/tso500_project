#!/usr/bin/env python3
"""
Combined OncoKB annotator for VEP-annotated VCF files (v3.0 — SMART).

- SNVs / indels / MantaINS  → /annotate/mutations/byProteinChange
- CNVs (MantaDUP/DEL)       → /annotate/copyNumberAlterations

Changes from v2.5:
  - CSQ field indices are now parsed dynamically from the VCF header,
    making the script compatible with any VEP version.
  - Uses shared csq_parser module (or inline fallback).
  - _extract_one_letter_protein() accepts hgvsp_index parameter.
  - find_preferred_csq_and_protein() accepts CSQ index parameters.
  - Guards against missing CSQ fields (index = -1).

Requires:
  - OncoKB token provided with --oncokb_token
  - VCF annotated by VEP with CSQ field
  - Optional preferred transcript list file for transcript prioritization
    (matches NM IDs ignoring version numbers)

Notes:
  - Selected OncoKB fields are flattened into dedicated INFO tags
  - The full raw OncoKB JSON is also stored in ONCOKB_JSON
"""

import os
import re
import json
import sys
import argparse
import requests
from cyvcf2 import VCF, Writer


# ---------------------------------------------------------
# Dynamic CSQ header parsing
# ---------------------------------------------------------

def parse_csq_format(vcf_path: str) -> tuple[list[str], dict[str, int]]:
    """
    Read a VCF file header and extract the CSQ field names and their indices.

    Returns:
        (field_names, field_map)
        - field_names: list of CSQ field names in order
        - field_map:   dict mapping field name → positional index
    """
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format:\s*(.+?)"', line)
                if m:
                    fields = m.group(1).split("|")
                    field_map = {name: i for i, name in enumerate(fields)}
                    return fields, field_map
            if line.startswith("#CHROM"):
                break

    raise RuntimeError(
        f"CSQ header (##INFO=<ID=CSQ,...>) not found in {vcf_path}. "
        "Is the VCF annotated by VEP?"
    )


def get_csq_index(field_map: dict[str, int], field_name: str, required: bool = False) -> int:
    """
    Look up a CSQ field index by name. Returns -1 if not found and not required.
    """
    if field_name in field_map:
        return field_map[field_name]
    if required:
        raise KeyError(
            f"Required CSQ field '{field_name}' not found in VEP annotation. "
            f"Available fields: {', '.join(sorted(field_map.keys()))}"
        )
    return -1


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

def load_preferred_transcripts(file_path: str | None) -> set:
    """
    Load a set of preferred NM IDs from a file, stripping version numbers.
    If no file is provided, return an empty set and print a warning.
    If the file does not exist or cannot be read, return an empty set and print a warning.
    This allows the script to continue and fall back to the default transcript selection logic.
    """
    if not file_path:
        sys.stderr.write(
            "WARNING: No preferred transcripts file was provided. "
            "Proceeding without transcript prioritization.\n"
        )
        return set()

    if not os.path.exists(file_path):
        sys.stderr.write(
            f"WARNING: Preferred transcripts file not found at {file_path}. "
            "Proceeding without transcript prioritization.\n"
        )
        return set()

    try:
        with open(file_path, "r") as f:
            return {
                line.strip().split(".")[0]
                for line in f
                if line.strip() and not line.startswith("#")
            }
    except Exception as e:
        sys.stderr.write(
            f"WARNING: Could not read preferred transcripts file '{file_path}': {e}. "
            "Proceeding without transcript prioritization.\n"
        )
        return set()


def clean_value(v):
    if v is None:
        return ""
    s = str(v)
    return s.replace("\n", " ").replace("\r", " ").replace(";", ",").strip()


def stringify_list(values) -> str:
    if not values:
        return ""
    out = []
    for v in values:
        cv = clean_value(v)
        if cv != "":
            out.append(cv)
    return "|".join(out)


def get_cancer_from_filename(filename: str) -> str:
    base = os.path.basename(filename)
    for part in re.split(r"[\_.\-]+", base):
        key = part.lower()
        if key in CANCER_MAP:
            return CANCER_MAP[key]
    return "UNKNOWN"


def _extract_one_letter_protein(csq: str, hgvsp_index: int):
    """
    Extract HGVSp from a CSQ annotation string and convert 3-letter amino
    acid codes to 1-letter codes where possible.
    """
    if hgvsp_index < 0:
        return None

    fields = csq.split("|")
    try:
        protein = fields[hgvsp_index]
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

    m = re.match(r"([A-Za-z\*]{3})(\d+)([A-Za-z\*]{3})", protein)
    if m:
        ref3, pos, alt3 = m.groups()
        ref1 = aa3to1.get(ref3, ref3[0])
        alt1 = aa3to1.get(alt3, alt3[0])
        return f"{ref1}{pos}{alt1}"

    return protein


def find_preferred_csq_and_protein(
    variant,
    csqs: str,
    preferred_tx: set,
    symbol_idx: int,
    feature_idx: int,
    hgvsp_idx: int,
    mane_select_nm_idx: int,
    mane_plus_clinical_nm_idx: int,
    mane_status_idx: int,
):
    """
    Finds the best CSQ annotation using a 3-tier prioritization logic.

    Returns: (gene, protein_change, selected_transcript, selection_reason)
    """
    if not csqs:
        return None, None, None, "No CSQ provided"

    annotations = csqs.split(",")
    selected_annotation = None
    selection_reason = "Fallback (VEP Default)"

    # --- TIER 1: PREFERRED TRANSCRIPTS (Match against NM IDs ignoring version) ---
    if preferred_tx:
        for ann in annotations:
            fields = ann.split("|")

            # Check MANE_SELECT NM ID
            if mane_select_nm_idx >= 0 and len(fields) > mane_select_nm_idx:
                current_nm_id_select = fields[mane_select_nm_idx].split(".")[0].strip()
                if current_nm_id_select and current_nm_id_select in preferred_tx:
                    selected_annotation = ann
                    selection_reason = "Preferred List Match (NM/MANE Select)"
                    break

            # Check MANE_PLUS_CLINICAL NM ID
            if (
                selected_annotation is None
                and mane_plus_clinical_nm_idx >= 0
                and len(fields) > mane_plus_clinical_nm_idx
            ):
                current_nm_id_plus = fields[mane_plus_clinical_nm_idx].split(".")[0].strip()
                if current_nm_id_plus and current_nm_id_plus in preferred_tx:
                    selected_annotation = ann
                    selection_reason = "Preferred List Match (NM/MANE Plus Clinical)"
                    break

    # --- TIER 2: MANE TRANSCRIPTS ---
    if selected_annotation is None and mane_status_idx >= 0:
        for ann in annotations:
            fields = ann.split("|")
            if len(fields) > mane_status_idx and fields[mane_status_idx] == "MANE Select":
                selected_annotation = ann
                selection_reason = "MANE Select"
                break

        if selected_annotation is None:
            for ann in annotations:
                fields = ann.split("|")
                if (
                    len(fields) > mane_status_idx
                    and fields[mane_status_idx].startswith("MANE Plus Clinical")
                ):
                    selected_annotation = ann
                    selection_reason = "MANE Plus Clinical"
                    break

    # --- TIER 3: FALLBACK ---
    if selected_annotation is None and annotations:
        selected_annotation = annotations[0]
        selection_reason = "Fallback (VEP Default)"

    if selected_annotation is None:
        selection_reason = "No valid CSQ found"
        return None, None, None, selection_reason

    fields = selected_annotation.split("|")
    gene = fields[symbol_idx] if len(fields) > symbol_idx else None
    selected_transcript = fields[feature_idx] if len(fields) > feature_idx else "N/A"

    protein_change = _extract_one_letter_protein(selected_annotation, hgvsp_idx)

    alt_string = ",".join(map(str, variant.ALT)) if variant.ALT else "."
    sys.stdout.write(
        f"  [Variant {variant.CHROM}:{variant.POS} {variant.REF}>{alt_string}] "
        f"Selected: {selection_reason}. Gene: {gene}. Tx: {selected_transcript}. Protein: {protein_change}\n"
    )

    return gene, protein_change, selected_transcript, selection_reason


# ---------------------------------------------------------
# OncoKB queries
# ---------------------------------------------------------

def query_oncokb_mutation(gene: str, alteration: str, tumor_type: str, token: str, cache: dict):
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


def query_oncokb_cna(gene: str, cna_type: str, tumor_type: str, token: str, cache: dict):
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

    query = data.get("query") or {}
    out["query_reference_genome"] = clean_value(query.get("referenceGenome"))
    out["query_hugo_symbol"] = clean_value(query.get("hugoSymbol"))
    out["query_entrez_gene_id"] = clean_value(query.get("entrezGeneId"))
    out["query_alteration"] = clean_value(query.get("alteration"))
    out["query_tumor_type"] = clean_value(query.get("tumorType"))

    out["gene_exist"] = clean_value(data.get("geneExist"))
    out["variant_exist"] = clean_value(data.get("variantExist"))
    out["allele_exist"] = clean_value(data.get("alleleExist"))
    out["oncogenic"] = clean_value(data.get("oncogenic"))
    out["hotspot"] = clean_value(data.get("hotspot"))
    out["vus"] = clean_value(data.get("VUS"))
    out["exon"] = clean_value(data.get("exon"))

    me = data.get("mutationEffect") or {}
    out["effect_known"] = clean_value(me.get("knownEffect"))
    out["effect_desc"] = clean_value(me.get("description"))
    citations = me.get("citations") or {}
    out["pmids"] = stringify_list(citations.get("pmids") or [])
    out["abstracts"] = stringify_list(citations.get("abstracts") or [])

    out["sens_level"] = clean_value(data.get("highestSensitiveLevel"))
    out["resist_level"] = clean_value(data.get("highestResistanceLevel"))
    out["diag_level"] = clean_value(data.get("highestDiagnosticImplicationLevel"))
    out["prog_level"] = clean_value(data.get("highestPrognosticImplicationLevel"))
    out["fda_level"] = clean_value(data.get("highestFdaLevel"))
    out["other_sensitive_levels"] = stringify_list(data.get("otherSignificantSensitiveLevels") or [])
    out["other_resistance_levels"] = stringify_list(data.get("otherSignificantResistanceLevels") or [])

    out["gene_summary"] = clean_value(data.get("geneSummary"))
    out["variant_summary"] = clean_value(data.get("variantSummary"))
    out["tumor_type_summary"] = clean_value(data.get("tumorTypeSummary"))
    out["diagnostic_summary"] = clean_value(data.get("diagnosticSummary"))
    out["prognostic_summary"] = clean_value(data.get("prognosticSummary"))

    diagnostic_implications = data.get("diagnosticImplications") or []
    prognostic_implications = data.get("prognosticImplications") or []
    treatments = data.get("treatments") or []
    out["diagnostic_implications"] = clean_value(json.dumps(diagnostic_implications, separators=(",", ":")))
    out["prognostic_implications"] = clean_value(json.dumps(prognostic_implications, separators=(",", ":")))
    out["treatments"] = clean_value(json.dumps(treatments, separators=(",", ":")))

    out["dataVersion"] = clean_value(data.get("dataVersion"))
    out["lastUpdate"] = clean_value(data.get("lastUpdate"))
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
      "snv": For SNVs, INDELs, MantaINS, MantaBND (treated as mutations).
    """
    svtype = variant.INFO.get("SVTYPE")
    variant_id = variant.ID if variant.ID else ""
    sv_class = "snv"

    if svtype:
        sv = str(svtype).upper()
        if sv in {"DUP", "DEL", "GAIN", "LOSS"}:
            sv_class = "cnv"
        elif sv in {"INS", "BND", "INV"}:
            sv_class = "snv"

    if "MantaDUP" in variant_id or "MantaDEL" in variant_id:
        sv_class = "cnv"
    elif "MantaINS" in variant_id or "MantaBND" in variant_id:
        sv_class = "snv"

    return sv_class


def map_svtype_to_cna(svtype: str, variant_id: str) -> str | None:
    """Map SVTYPE or ID → OncoKB copyNameAlterationType."""
    if variant_id:
        if "MantaDUP" in variant_id:
            return "AMPLIFICATION"
        if "MantaDEL" in variant_id:
            return "DELETION"

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

def annotate_vcf(
    input_vcf: str,
    output_vcf: str,
    tumor_mode: str,
    oncokb_token: str,
    preferred_transcripts_file: str | None = None,
):
    token = oncokb_token
    preferred_transcripts = load_preferred_transcripts(preferred_transcripts_file)

    if tumor_mode == "filename":
        tumor_type = get_cancer_from_filename(input_vcf)
    elif tumor_mode == "generic":
        tumor_type = "UNKNOWN"
    else:
        sys.stderr.write(f"ERROR: Invalid tumor mode '{tumor_mode}'.\n")
        sys.exit(1)

    if tumor_type == "UNKNOWN" and tumor_mode == "filename":
        sys.stderr.write(
            f"WARNING: Tumor type could not be determined from filename. Using '{tumor_type}'.\n"
        )

    print(f"\n✔ Annotation Mode: {'Filename Inference' if tumor_mode == 'filename' else 'Generic Query'}")
    print(f"✔ Detected tumor type used for OncoKB: {tumor_type}")

    if preferred_transcripts:
        print(f"✔ Using {len(preferred_transcripts)} preferred NM transcripts (version agnostic).")
    else:
        print("⚠️  No preferred transcripts loaded. Falling back to MANE/VEP default transcript selection.")

    # --- Dynamically parse CSQ field positions from the VCF header ---
    _csq_fields, _csq_map = parse_csq_format(input_vcf)

    SYMBOL_IDX                = get_csq_index(_csq_map, "SYMBOL", required=True)
    FEATURE_IDX               = get_csq_index(_csq_map, "Feature", required=True)
    HGVSP_IDX                 = get_csq_index(_csq_map, "HGVSp", required=True)
    MANE_SELECT_NM_IDX        = get_csq_index(_csq_map, "MANE_SELECT")
    MANE_PLUS_CLINICAL_NM_IDX = get_csq_index(_csq_map, "MANE_PLUS_CLINICAL")
    MANE_STATUS_IDX           = get_csq_index(_csq_map, "MANE")

    print(f"✔ Parsed {len(_csq_fields)} CSQ fields from VCF header")
    print(
        f"  SYMBOL={SYMBOL_IDX}, Feature={FEATURE_IDX}, HGVSp={HGVSP_IDX}, "
        f"MANE_SELECT={MANE_SELECT_NM_IDX}, MANE_PLUS_CLINICAL={MANE_PLUS_CLINICAL_NM_IDX}, "
        f"MANE={MANE_STATUS_IDX}"
    )

    vcf = VCF(input_vcf)

    fields_to_add = {
        "ONCOKB_QUERY_REF_GENOME": "OncoKB query.referenceGenome",
        "ONCOKB_QUERY_HUGO_SYMBOL": "OncoKB query.hugoSymbol",
        "ONCOKB_QUERY_ENTREZ_GENE_ID": "OncoKB query.entrezGeneId",
        "ONCOKB_QUERY_ALTERATION": "OncoKB query.alteration",
        "ONCOKB_QUERY_TUMOR_TYPE": "OncoKB query.tumorType",
        "ONCOKB_ONCOGENIC": "OncoKB oncogenic classification",
        "ONCOKB_HOTSPOT": "OncoKB hotspot flag",
        "ONCOKB_EFFECT": "OncoKB mutationEffect.knownEffect",
        "ONCOKB_EFFECT_DESC": "OncoKB mutationEffect.description",
        "ONCOKB_SENS_LVL": "OncoKB highestSensitiveLevel",
        "ONCOKB_RESIST_LVL": "OncoKB highestResistanceLevel",
        "ONCOKB_DIAG_LVL": "OncoKB highestDiagnosticImplicationLevel",
        "ONCOKB_PROG_LVL": "OncoKB highestPrognosticImplicationLevel",
        "ONCOKB_FDA_LVL": "OncoKB highestFdaLevel",
        "ONCOKB_OTHER_SENS_LVLS": "OncoKB otherSignificantSensitiveLevels",
        "ONCOKB_OTHER_RESIST_LVLS": "OncoKB otherSignificantResistanceLevels",
        "ONCOKB_TREATMENTS": "OncoKB treatments JSON",
        "ONCOKB_GENE_SUMMARY": "OncoKB geneSummary",
        "ONCOKB_VARIANT_SUMMARY": "OncoKB variantSummary",
        "ONCOKB_TUMOR_TYPE_SUMMARY": "OncoKB tumorTypeSummary",
        "ONCOKB_DIAG_SUMMARY": "OncoKB diagnosticSummary",
        "ONCOKB_PROG_SUMMARY": "OncoKB prognosticSummary",
        "ONCOKB_VUS": "OncoKB VUS classification",
        "ONCOKB_GENE_EXIST": "OncoKB geneExist",
        "ONCOKB_VARIANT_EXIST": "OncoKB variantExist",
        "ONCOKB_ALLELE_EXIST": "OncoKB alleleExist",
        "ONCOKB_EXON": "OncoKB exon",
        "ONCOKB_PMIDS": "OncoKB mutationEffect citations pmids",
        "ONCOKB_ABSTRACTS": "OncoKB mutationEffect citations abstracts",
        "ONCOKB_DIAGNOSTIC_IMPLICATIONS": "OncoKB diagnosticImplications JSON",
        "ONCOKB_PROGNOSTIC_IMPLICATIONS": "OncoKB prognosticImplications JSON",
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
        v_class = classify_variant(variant)
        csqs = variant.INFO.get("CSQ")
        gene = None
        protein_change = None
        summary = {}

        if v_class == "snv":
            gene, protein_change, selected_transcript, selection_reason = find_preferred_csq_and_protein(
                variant, str(csqs), preferred_transcripts,
                symbol_idx=SYMBOL_IDX,
                feature_idx=FEATURE_IDX,
                hgvsp_idx=HGVSP_IDX,
                mane_select_nm_idx=MANE_SELECT_NM_IDX,
                mane_plus_clinical_nm_idx=MANE_PLUS_CLINICAL_NM_IDX,
                mane_status_idx=MANE_STATUS_IDX,
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
            if csqs:
                fields = str(csqs).split(",")[0].split("|")
                gene = fields[SYMBOL_IDX] if len(fields) > SYMBOL_IDX else None

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
            "query_reference_genome": "ONCOKB_QUERY_REF_GENOME",
            "query_hugo_symbol": "ONCOKB_QUERY_HUGO_SYMBOL",
            "query_entrez_gene_id": "ONCOKB_QUERY_ENTREZ_GENE_ID",
            "query_alteration": "ONCOKB_QUERY_ALTERATION",
            "query_tumor_type": "ONCOKB_QUERY_TUMOR_TYPE",
            "oncogenic": "ONCOKB_ONCOGENIC",
            "hotspot": "ONCOKB_HOTSPOT",
            "vus": "ONCOKB_VUS",
            "gene_exist": "ONCOKB_GENE_EXIST",
            "variant_exist": "ONCOKB_VARIANT_EXIST",
            "allele_exist": "ONCOKB_ALLELE_EXIST",
            "exon": "ONCOKB_EXON",
            "dataVersion": "ONCOKB_DATA_VERSION",
            "lastUpdate": "ONCOKB_LAST_UPDATE",
            "effect_known": "ONCOKB_EFFECT",
            "effect_desc": "ONCOKB_EFFECT_DESC",
            "pmids": "ONCOKB_PMIDS",
            "abstracts": "ONCOKB_ABSTRACTS",
            "sens_level": "ONCOKB_SENS_LVL",
            "resist_level": "ONCOKB_RESIST_LVL",
            "diag_level": "ONCOKB_DIAG_LVL",
            "prog_level": "ONCOKB_PROG_LVL",
            "fda_level": "ONCOKB_FDA_LVL",
            "other_sensitive_levels": "ONCOKB_OTHER_SENS_LVLS",
            "other_resistance_levels": "ONCOKB_OTHER_RESIST_LVLS",
            "gene_summary": "ONCOKB_GENE_SUMMARY",
            "variant_summary": "ONCOKB_VARIANT_SUMMARY",
            "tumor_type_summary": "ONCOKB_TUMOR_TYPE_SUMMARY",
            "diagnostic_summary": "ONCOKB_DIAG_SUMMARY",
            "prognostic_summary": "ONCOKB_PROG_SUMMARY",
            "diagnostic_implications": "ONCOKB_DIAGNOSTIC_IMPLICATIONS",
            "prognostic_implications": "ONCOKB_PROGNOSTIC_IMPLICATIONS",
            "treatments": "ONCOKB_TREATMENTS",
            "json": "ONCOKB_JSON",
        }

        for key, value in summary.items():
            if key in mapping and value != "":
                variant.INFO[mapping[key]] = value

        writer.write_record(variant)

    writer.close()

    print("\n📊 OncoKB Annotation Summary:")
    print(f"  Total variants processed: {total_variants}")
    print(f"  Mutations annotated:      {annotated_mut}")
    print(f"  CNVs annotated:           {annotated_cna}")
    if total_variants:
        print(f"  Overall annotation rate:  {(annotated_mut + annotated_cna) / total_variants * 100:.2f}%")
    print(f"\n✔ OncoKB-annotated VCF written to: {output_vcf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate VCF files using OncoKB.")
    parser.add_argument("input_vcf", help="Input VCF file annotated by VEP.")
    parser.add_argument("output_vcf", help="Output VCF file.")
    parser.add_argument(
        "--tumor_mode",
        choices=["filename", "generic"],
        default="filename",
        help="Tumor type logic: 'filename' (infer from name) or 'generic' (UNKNOWN).",
    )
    parser.add_argument(
        "--oncokb_token",
        required=True,
        help="Token to get access to the OncoKB API.",
    )
    parser.add_argument(
        "--preferred-transcripts",
        default=None,
        help=(
            "Optional path to a preferred transcript list file. "
            "If not provided, the script will warn and continue without transcript prioritization."
        ),
    )
    args = parser.parse_args()

    annotate_vcf(
        args.input_vcf,
        args.output_vcf,
        args.tumor_mode,
        args.oncokb_token,
        args.preferred_transcripts,
    )
