#!/usr/bin/env python3
"""
SMART pipeline verification — VEP + OncoKB field-level comparison.

Runs both verifications against a Final_result_tier1.maf and reports PASS /
MISMATCH for every checked field.  Results are written to a single TSV with a
'Module' column ('VEP' or 'OncoKB') so downstream analysis is easy.

Usage:
    # Run both modules (default)
    python tests/verification_combined/verify.py \
        --maf  tests/verification_combined/output/output/Final_result_tier1.maf \
        --token $ONCOKB_TOKEN \
        --output results.tsv

    # Run only OncoKB
    python tests/verification_combined/verify.py \
        --maf tests/verification_combined/output/output/Final_result_tier1.maf \
        --token $ONCOKB_TOKEN \
        --modules oncokb

    # Run only VEP (no token needed)
    python tests/verification_combined/verify.py \
        --maf tests/verification_combined/output/output/Final_result_tier1.maf \
        --modules vep

    # Limit to specific variant IDs
    python tests/verification_combined/verify.py \
        --maf ... --token ... \
        --ids BRAF_V600E IDH1_R132H MantaDUP:ERBB2_AMP

Exit code 0 if all checks pass, 1 if any MISMATCH or ERROR.

Fields NOT checked (require local databases/plugins not available via REST API):
  gnomAD*, ClinVar*, SpliceAI*, REVEL, LOEUF, CancerHotspots*, CIViC*
"""

import argparse
import re
import sys
import time
import urllib.parse
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

# =============================================================================
# Shared utilities
# =============================================================================

def load_maf(path: Path) -> pd.DataFrame:
    with open(path) as fh:
        first = fh.readline()
    skip = 1 if first.startswith("#") else 0
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str,
                       low_memory=False).fillna("")


def _norm(v) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s in ("nan", "None", "NA", "N/A") else s


def maf_value(row: pd.Series, col: str) -> str:
    return _norm(row.get(col, ""))


def is_cna(row: pd.Series) -> bool:
    vid = str(row.get("ID", ""))
    svtype = str(row.get("SVTYPE", ""))
    return (vid.startswith("MantaDUP:") or vid.startswith("MantaDEL:")
            or svtype in ("DUP", "DEL"))


def print_coverage_report(results: list, field_map: list,
                           skip_sentinel_fields: set,
                           hints: dict, field_map_width: int = 3) -> None:
    """
    Print a coverage summary.
    field_map_width: number of elements per FIELD_MAP tuple (3 for OncoKB, 4 for VEP).
    """
    field_data = defaultdict(list)
    for r in results:
        if r.get("MAF_Field") in skip_sentinel_fields:
            continue
        field_data[r["MAF_Field"]].append((r["MAF_Value"], r["API_Value"]))

    untestable, tested = [], []
    for entry in field_map:
        maf_col = entry[0]
        pairs = field_data.get(maf_col, [])
        if any(mv or av for mv, av in pairs):
            tested.append(maf_col)
        else:
            untestable.append(maf_col)

    print(f"\n{'─'*60}")
    print("FIELD COVERAGE REPORT")
    print(f"{'─'*60}")
    print(f"  Tested (≥1 variant had data):  {len(tested)}")
    for f in tested:
        print(f"    ✓ {f}")
    if untestable:
        print(f"\n  UNTESTABLE (always empty in MAF + API): {len(untestable)}")
        for f in untestable:
            hint = hints.get(f, "No hint available.")
            print(f"    ○ {f}")
            print(f"        → {hint}")
    print(f"{'─'*60}")


# =============================================================================
# OncoKB verification
# =============================================================================

ONCOKB_BASE = "https://www.oncokb.org/api/v1"
REF_GENOME  = "GRCh38"

# --- helpers ---

def _str(v):
    return "" if v is None else str(v).strip()

def _bool_str(v):
    if v is True:  return "True"
    if v is False: return "False"
    return ""

def _treatments_by_level(data, level):
    drugs = []
    for t in (data.get("treatments") or []):
        if t.get("level") == level:
            for d in t.get("drugs", []):
                name = d.get("drugName", "")
                if name and name not in drugs:
                    drugs.append(name)
    return ",".join(drugs)

def _normalise_drugs(val: str) -> set:
    if not val:
        return set()
    drugs = set()
    for part in val.replace("+", ",").split(","):
        name = part.strip()
        if name:
            drugs.add(name.lower())
    return drugs

DRUG_FIELDS = {
    "LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "LEVEL_R1", "LEVEL_R2",
}

# (maf_col, api_extractor, label)
ONCOKB_FIELD_MAP = [
    # NOTE: ONCOKB_LAST_UPDATE / ONCOKB_DATA_VERSION excluded — MafAnnotator
    # (batch endpoint) and byProteinChange (single endpoint) return slightly
    # different timestamps for the same variant. API endpoint difference, not a
    # pipeline bug.
    ("ONCOKB_ALLELE_EXIST",     lambda d: _bool_str(d.get("alleleExist")),                         "alleleExist"),
    ("ONCOKB_GENE_SUMMARY",     lambda d: _str(d.get("geneSummary")),                              "geneSummary"),
    ("ONCOKB_HOTSPOT",          lambda d: _bool_str(d.get("hotspot")),                             "hotspot"),
    ("ONCOKB_VUS",              lambda d: _bool_str(d.get("VUS")),                                 "VUS"),
    ("ONCOKB_highestFdaLevel",  lambda d: _str(d.get("highestFdaLevel") or ""),                    "highestFdaLevel"),
    ("ONCOKB_VARIANT_SUMMARY",  lambda d: _str(d.get("variantSummary")),                           "variantSummary"),
    ("GENE_IN_ONCOKB",          lambda d: _bool_str(d.get("geneExist")),                           "geneExist"),
    ("VARIANT_IN_ONCOKB",       lambda d: _bool_str(d.get("variantExist")),                        "variantExist"),
    ("MUTATION_EFFECT",         lambda d: _str((d.get("mutationEffect") or {}).get("knownEffect")), "mutationEffect.knownEffect"),
    # MUTATION_EFFECT_CITATIONS: MafAnnotator includes abstracts + extra PMIDs
    # beyond what the single endpoint returns — subset only. Skipped.
    ("ONCOGENIC",               lambda d: _str(d.get("oncogenic")),                                "oncogenic"),
    ("LEVEL_1",                 lambda d: _treatments_by_level(d, "LEVEL_1"),                      "treatments LEVEL_1"),
    ("LEVEL_2",                 lambda d: _treatments_by_level(d, "LEVEL_2"),                      "treatments LEVEL_2"),
    ("LEVEL_3A",                lambda d: _treatments_by_level(d, "LEVEL_3A"),                     "treatments LEVEL_3A"),
    ("LEVEL_3B",                lambda d: _treatments_by_level(d, "LEVEL_3B"),                     "treatments LEVEL_3B"),
    ("LEVEL_4",                 lambda d: _treatments_by_level(d, "LEVEL_4"),                      "treatments LEVEL_4"),
    ("LEVEL_R1",                lambda d: _treatments_by_level(d, "LEVEL_R1"),                     "treatments LEVEL_R1"),
    ("LEVEL_R2",                lambda d: _treatments_by_level(d, "LEVEL_R2"),                     "treatments LEVEL_R2"),
    # HIGHEST_LEVEL includes resistance levels in ranking (R1 > L4) — not
    # directly comparable to any single API field. Skipped.
    ("HIGHEST_SENSITIVE_LEVEL",  lambda d: _str(d.get("highestSensitiveLevel")),                   "highestSensitiveLevel"),
    ("HIGHEST_RESISTANCE_LEVEL", lambda d: _str(d.get("highestResistanceLevel")),                  "highestResistanceLevel"),
    ("ONCOKB_DIAG_LVL",          lambda d: _str(d.get("highestDiagnosticImplicationLevel") or ""), "highestDiagnosticImplicationLevel"),
    ("ONCOKB_PROG_LVL",          lambda d: _str(d.get("highestPrognosticImplicationLevel") or ""), "highestPrognosticImplicationLevel"),
    ("ONCOKB_diagnosticSummary", lambda d: _str(d.get("diagnosticSummary")),                       "diagnosticSummary"),
    ("ONCOKB_prognosticSummary", lambda d: _str(d.get("prognosticSummary")),                       "prognosticSummary"),
]

ONCOKB_COVERAGE_HINTS = {
    "ONCOKB_diagnosticSummary": (
        "Always empty when querying without a tumorType. The API returns "
        "highestDiagnosticImplicationLevel (ONCOKB_DIAG_LVL) pan-cancer, but "
        "diagnosticSummary text is only populated for a specific tumor type. "
        "Variants in this dataset that HAVE a DIAG level: IDH1 R132H (Dx2), "
        "BRAF V600E (Dx2), NRAS Q61R (Dx2), KRAS G12D/G13D (Dx2), CDKN2A_DEL (Dx2), "
        "PTEN R130* (Dx3). To test, re-query with e.g. tumorType=Glioma for IDH1 R132H."
    ),
    "ONCOKB_prognosticSummary": (
        "Always empty when querying without a tumorType. Same reason as "
        "diagnosticSummary. Variants with PROG level: IDH1 R132H (Px2), "
        "NRAS Q61R (Px1), TP53 R175H (Px1). To test, re-query with a specific tumor type."
    ),
}


def clean_protein_change(hgvsp_short: str) -> str:
    s = hgvsp_short.strip()
    if ":" in s:
        s = s.split(":")[-1]
    if s.startswith("p."):
        s = s[2:]
    return s


def query_oncokb_mutation(hugo: str, alteration: str, token: str) -> dict:
    url = f"{ONCOKB_BASE}/annotate/mutations/byProteinChange"
    params = {"hugoSymbol": hugo, "alteration": alteration,
               "referenceGenome": REF_GENOME}
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}
    r = requests.get(url, params=params, headers=headers, timeout=15)
    r.raise_for_status()
    return r.json()


def query_oncokb_cna(hugo: str, cna_type: str, token: str) -> dict:
    url = f"{ONCOKB_BASE}/annotate/copyNumberAlterations"
    params = {"hugoSymbol": hugo, "copyNameAlterationType": cna_type,
               "referenceGenome": REF_GENOME}
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}
    r = requests.get(url, params=params, headers=headers, timeout=15)
    r.raise_for_status()
    return r.json()


def compare_oncokb(maf_val: str, api_val: str, maf_col: str) -> str:
    if maf_val == api_val:                        return "PASS"
    if maf_val.lower() == api_val.lower():        return "PASS"
    if maf_col in DRUG_FIELDS:
        if _normalise_drugs(maf_val) == _normalise_drugs(api_val):
            return "PASS"
    return "MISMATCH"


def run_oncokb(df: pd.DataFrame, token: str) -> list:
    results = []
    total_pass = total_mismatch = total_error = 0

    print(f"\n{'#'*60}")
    print("MODULE: OncoKB")
    print(f"{'#'*60}")

    for i, (_, row) in enumerate(df.iterrows()):
        hugo   = str(row.get("Hugo_Symbol", "")).strip()
        vid    = str(row.get("ID", f"row_{i}")).strip()
        hgvsp  = str(row.get("HGVSp_Short", "")).strip()
        sample = str(row.get("Tumor_Sample_Barcode", "")).strip()

        print(f"\n{'─'*60}")
        print(f"[{i+1}/{len(df)}]  {sample}  |  {hugo}  |  {vid}  |  {hgvsp}")

        try:
            if is_cna(row):
                vid_upper = vid.upper()
                cna_type = "AMPLIFICATION" if ("DUP" in vid_upper or "AMP" in vid_upper) else "DELETION"
                print(f"  → CNA query: {hugo} {cna_type}")
                api_data = query_oncokb_cna(hugo, cna_type, token)
            else:
                alteration = clean_protein_change(hgvsp)
                if not alteration:
                    print(f"  SKIP: no protein change for {vid}")
                    continue
                print(f"  → Mutation query: {hugo} {alteration}")
                api_data = query_oncokb_mutation(hugo, alteration, token)
            time.sleep(0.2)
        except Exception as exc:
            print(f"  ERROR querying OncoKB: {exc}")
            total_error += 1
            results.append({
                "Module": "OncoKB", "Sample": sample, "Variant_ID": vid,
                "Gene": hugo, "HGVSp": hgvsp, "MAF_Field": "API_CALL",
                "API_Field": "", "MAF_Value": "", "API_Value": "",
                "Result": f"ERROR: {exc}", "Note": "",
            })
            continue

        for maf_col, api_fn, label in ONCOKB_FIELD_MAP:
            maf_val = maf_value(row, maf_col)
            try:
                api_val = _norm(api_fn(api_data))
            except Exception as exc:
                api_val = f"EXTRACT_ERROR:{exc}"

            status = compare_oncokb(maf_val, api_val, maf_col)
            icon   = "✓" if status == "PASS" else "✗"
            print(f"  {icon} {maf_col:<35} MAF={maf_val[:55]!r}  API={api_val[:55]!r}  [{status}]")

            if status == "PASS": total_pass += 1
            else:                total_mismatch += 1

            results.append({
                "Module":      "OncoKB",
                "Sample":      sample,
                "Variant_ID":  vid,
                "Gene":        hugo,
                "HGVSp":       hgvsp,
                "MAF_Field":   maf_col,
                "API_Field":   label,
                "MAF_Value":   maf_val,
                "API_Value":   api_val,
                "Result":      status,
                "Note":        "",
            })

    # Summary
    total = total_pass + total_mismatch + total_error
    print(f"\n{'='*60}")
    print("OncoKB SUMMARY")
    print(f"{'='*60}")
    print(f"  Variants checked:  {len(df)}")
    print(f"  Field checks:      {total}")
    print(f"  PASS:              {total_pass}")
    print(f"  MISMATCH:          {total_mismatch}")
    print(f"  ERROR:             {total_error}")
    print(f"{'='*60}")

    mismatches = [r for r in results if r["Result"] not in ("PASS",)]
    if mismatches:
        print("\nMISMATCHES / ERRORS:")
        for m in mismatches:
            print(f"  [{m['Variant_ID']}] {m['MAF_Field']}")
            print(f"      MAF: {m['MAF_Value'][:80]!r}")
            print(f"      API: {m['API_Value'][:80]!r}")

    print_coverage_report(results, ONCOKB_FIELD_MAP,
                          skip_sentinel_fields={"API_CALL"},
                          hints=ONCOKB_COVERAGE_HINTS,
                          field_map_width=3)
    return results


# =============================================================================
# VEP verification
# =============================================================================

ENSEMBL_BASE = "https://rest.ensembl.org"
VEP_HEADERS  = {"Content-Type": "application/json", "Accept": "application/json"}


def _strip_version(tid: str) -> str:
    if tid.startswith(("NM_", "NR_", "NP_")):
        return tid
    return tid.split(".")[0]

def _strip_hgvsc_prefix(hgvsc: str) -> str:
    return hgvsc.split(":")[-1] if ":" in hgvsc else hgvsc

def _strip_hgvsp_version(hgvsp: str) -> str:
    if ":" in hgvsp:
        prot_id, notation = hgvsp.split(":", 1)
        return f"{prot_id.split('.')[0]}:{notation}"
    return hgvsp

def _consequences(t: dict) -> str:
    return ",".join(sorted(t.get("consequence_terms", [])))

def _sift_str(t: dict) -> str:
    pred  = t.get("sift_prediction", "")
    score = t.get("sift_score")
    if not pred: return ""
    return f"{pred}({score})" if score is not None else pred

def _polyphen_str(t: dict) -> str:
    pred  = t.get("polyphen_prediction", "")
    score = t.get("polyphen_score")
    if not pred: return ""
    return f"{pred}({score})" if score is not None else pred

def _protein_pos(t: dict) -> str:
    start = t.get("protein_start")
    end   = t.get("protein_end")
    if start is None: return ""
    if end is not None and end != start: return f"{start}-{end}"
    return str(start)

# (maf_col, api_extractor, label, note)
VEP_FIELD_MAP = [
    ("Consequence",      _consequences,                                       "consequence_terms",           "sorted, comma-joined"),
    ("IMPACT",           lambda t: t.get("impact", ""),                      "impact",                      ""),
    ("BIOTYPE",          lambda t: t.get("biotype", ""),                     "biotype",                     ""),
    ("EXON",             lambda t: t.get("exon", ""),                        "exon",                        "e.g. 15/18"),
    ("INTRON",           lambda t: t.get("intron", ""),                      "intron",                      "only for intronic variants"),
    ("HGVSc",            lambda t: _strip_hgvsc_prefix(t.get("hgvsc", "")), "hgvsc (c. part only)",        "transcript version stripped"),
    ("HGVSp",            lambda t: _strip_hgvsp_version(t.get("hgvsp", "")),"hgvsp (version stripped)",    "3-letter AA notation"),
    ("Protein_position", _protein_pos,                                       "protein_start/end",            ""),
    ("Amino_acids",      lambda t: t.get("amino_acids", ""),                 "amino_acids",                  "e.g. V/E"),
    ("Codons",           lambda t: t.get("codons", ""),                      "codons",                       "case-insensitive"),
    ("SIFT",             _sift_str,                                           "sift_prediction(score)",       "score tolerance 1e-4"),
    ("PolyPhen",         _polyphen_str,                                       "polyphen_prediction(score)",   "score tolerance 1e-4"),
    ("STRAND",           lambda t: str(t.get("strand", "")),                 "strand",                       "1 or -1"),
    ("CANONICAL",        lambda t: "YES" if t.get("canonical") else "",      "canonical",                    ""),
]

VEP_COVERAGE_HINTS = {
    "INTRON": (
        "Only populated for intronic variants; no intronic variants in this dataset."
    ),
    "SIFT": (
        "Only populated for missense variants where SIFT has a model for the gene. "
        "NRAS Q61R, IDH1 R132H, BRAF V600E etc. should populate this field."
    ),
    "PolyPhen": (
        "Only populated for missense variants where PolyPhen has a model. "
        "Not all genes are covered."
    ),
}


def query_vep_hgvs(hgvs: str, refseq: bool = False) -> list:
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{ENSEMBL_BASE}/vep/human/hgvs/{encoded}"
    params = {
        "hgvs":      1,   # include HGVSc/HGVSp
        "canonical": 1,   # flag canonical transcript
        "sift":      "b", # prediction + score
        "polyphen":  "b", # prediction + score
        "numbers":   1,   # exon/intron numbers
    }
    if refseq:
        params["refseq"] = 1
    r = requests.get(url, params=params, headers=VEP_HEADERS, timeout=20)
    if r.status_code in (400, 404):
        return []
    r.raise_for_status()
    return r.json()


def find_transcript(vep_results: list, feature_id: str):
    feat_nv = _strip_version(feature_id)
    for result in vep_results:
        for t in result.get("transcript_consequences", []):
            tid = t.get("transcript_id", "")
            if tid == feature_id or _strip_version(tid) == feat_nv:
                return t
    for result in vep_results:
        tcs = result.get("transcript_consequences", [])
        if tcs:
            return tcs[0]
    return None


def _score_from(s: str):
    m = re.search(r'\(([0-9.e+\-]+)\)', s)
    if m:
        try: return float(m.group(1))
        except ValueError: pass
    return None


def compare_vep(maf_val: str, api_val: str, col: str) -> str:
    if maf_val == api_val:               return "PASS"
    if maf_val.lower() == api_val.lower(): return "PASS"
    if col == "Consequence":
        if set(maf_val.split(",")) == set(api_val.split(",")): return "PASS"
    if col == "HGVSp":
        if _strip_hgvsp_version(maf_val) == api_val: return "PASS"
    if col == "Protein_position":
        if maf_val.split("-")[0] == api_val.split("-")[0]: return "PASS"
    if col in ("SIFT", "PolyPhen"):
        ml = maf_val.split("(")[0].strip().lower()
        al = api_val.split("(")[0].strip().lower()
        if ml and ml == al:
            ms, as_ = _score_from(maf_val), _score_from(api_val)
            if ms is not None and as_ is not None:
                if abs(ms - as_) < 1e-4: return "PASS"
            elif ms is None and as_ is None:
                return "PASS"
    return "MISMATCH"


def run_vep(df: pd.DataFrame) -> list:
    results = []
    total_pass = total_mismatch = total_error = 0
    n_skipped = 0

    print(f"\n{'#'*60}")
    print("MODULE: VEP")
    print(f"{'#'*60}")

    for i, (_, row) in enumerate(df.iterrows()):
        hugo    = str(row.get("Hugo_Symbol", "")).strip()
        vid     = str(row.get("ID", f"row_{i}")).strip()
        feature = str(row.get("Feature", "")).strip()
        hgvsc   = str(row.get("HGVSc", "")).strip()
        hgvsp   = str(row.get("HGVSp_Short", "")).strip()
        sample  = str(row.get("Tumor_Sample_Barcode", "")).strip()

        print(f"\n{'─'*60}")
        print(f"[{i+1}/{len(df)}]  {sample}  |  {hugo}  |  {vid}  |  {hgvsp}")

        if is_cna(row):
            print("  SKIP: CNA — VEP REST API does not annotate structural variants")
            n_skipped += 1
            continue

        if not hgvsc or not feature:
            print(f"  SKIP: missing HGVSc ({hgvsc!r}) or Feature ({feature!r})")
            n_skipped += 1
            results.append({
                "Module": "VEP", "Sample": sample, "Variant_ID": vid,
                "Gene": hugo, "HGVSp": hgvsp, "MAF_Field": "NO_HGVSC",
                "API_Field": "", "MAF_Value": "", "API_Value": "",
                "Result": "SKIP: no HGVSc/Feature", "Note": "",
            })
            continue

        hgvs_query = f"{feature}:{hgvsc}"
        is_refseq  = feature.startswith(("NM_", "NR_"))
        print(f"  → VEP query: {hgvs_query}  (refseq={is_refseq})")

        try:
            vep_results = query_vep_hgvs(hgvs_query, refseq=is_refseq)
            time.sleep(0.5)
        except Exception as exc:
            print(f"  ERROR querying Ensembl VEP: {exc}")
            total_error += 1
            results.append({
                "Module": "VEP", "Sample": sample, "Variant_ID": vid,
                "Gene": hugo, "HGVSp": hgvsp, "MAF_Field": "API_CALL",
                "API_Field": "", "MAF_Value": "", "API_Value": "",
                "Result": f"ERROR: {exc}", "Note": "",
            })
            continue

        if not vep_results:
            print(f"  SKIP: VEP REST returned no results for {hgvs_query!r}")
            n_skipped += 1
            results.append({
                "Module": "VEP", "Sample": sample, "Variant_ID": vid,
                "Gene": hugo, "HGVSp": hgvsp, "MAF_Field": "API_CALL",
                "API_Field": "", "MAF_Value": "", "API_Value": "",
                "Result": "SKIP: no VEP results", "Note": "",
            })
            continue

        transcript = find_transcript(vep_results, feature)
        if transcript is None:
            print("  SKIP: no transcript consequences in VEP response")
            n_skipped += 1
            results.append({
                "Module": "VEP", "Sample": sample, "Variant_ID": vid,
                "Gene": hugo, "HGVSp": hgvsp, "MAF_Field": "NO_TRANSCRIPT",
                "API_Field": "", "MAF_Value": "", "API_Value": "",
                "Result": "SKIP: no transcript in VEP response", "Note": "",
            })
            continue

        used_tid = transcript.get("transcript_id", "?")
        if _strip_version(used_tid) != _strip_version(feature):
            print(f"  NOTE: exact transcript not found; using fallback {used_tid!r}")

        for maf_col, api_fn, label, note in VEP_FIELD_MAP:
            maf_val = maf_value(row, maf_col)
            try:
                api_val = _norm(api_fn(transcript))
            except Exception as exc:
                api_val = f"EXTRACT_ERROR:{exc}"

            status = compare_vep(maf_val, api_val, maf_col)
            icon   = "✓" if status == "PASS" else "✗"
            print(f"  {icon} {maf_col:<20} MAF={maf_val[:50]!r}  API={api_val[:50]!r}  [{status}]")

            if status == "PASS": total_pass += 1
            else:                total_mismatch += 1

            results.append({
                "Module":      "VEP",
                "Sample":      sample,
                "Variant_ID":  vid,
                "Gene":        hugo,
                "HGVSp":       hgvsp,
                "MAF_Field":   maf_col,
                "API_Field":   label,
                "MAF_Value":   maf_val,
                "API_Value":   api_val,
                "Result":      status,
                "Note":        note,
            })

    # Summary
    total = total_pass + total_mismatch + total_error
    print(f"\n{'='*60}")
    print("VEP SUMMARY")
    print(f"{'='*60}")
    print(f"  Rows in MAF:       {len(df)}")
    print(f"  Skipped (CNA/etc): {n_skipped}")
    print(f"  Variants checked:  {len(df) - n_skipped - total_error}")
    print(f"  Field checks:      {total}")
    print(f"  PASS:              {total_pass}")
    print(f"  MISMATCH:          {total_mismatch}")
    print(f"  ERROR:             {total_error}")
    print(f"{'='*60}")

    mismatches = [r for r in results if r["Result"] == "MISMATCH"]
    if mismatches:
        print("\nMISMATCHES:")
        for m in mismatches:
            print(f"  [{m['Variant_ID']}] {m['MAF_Field']}")
            print(f"      MAF: {m['MAF_Value'][:80]!r}")
            print(f"      API: {m['API_Value'][:80]!r}")

    print_coverage_report(results, VEP_FIELD_MAP,
                          skip_sentinel_fields={"API_CALL", "NO_HGVSC", "NO_TRANSCRIPT"},
                          hints=VEP_COVERAGE_HINTS,
                          field_map_width=4)
    return results


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Verify VEP and/or OncoKB fields in SMART MAF output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--maf",     required=True,
                        help="Path to Final_result_tier1.maf")
    parser.add_argument("--token",   default=None,
                        help="OncoKB API bearer token (required when running oncokb module)")
    parser.add_argument("--modules", default=["all"], nargs="+",
                        choices=["all", "vep", "oncokb"],
                        help="Which modules to run (default: all)")
    parser.add_argument("--output",  default=None,
                        help="Optional TSV output path (combined results from all modules)")
    parser.add_argument("--ids",     default=None, nargs="+",
                        help="Only check these variant IDs. Default: all rows.")
    args = parser.parse_args()

    modules = set(args.modules)
    if "all" in modules:
        modules = {"vep", "oncokb"}

    if "oncokb" in modules and not args.token:
        print("ERROR: --token is required when running the oncokb module", file=sys.stderr)
        sys.exit(1)

    maf_path = Path(args.maf)
    if not maf_path.exists():
        print(f"ERROR: MAF not found: {maf_path}", file=sys.stderr)
        sys.exit(1)

    df = load_maf(maf_path)
    print(f"Loaded MAF: {maf_path.name}  ({len(df)} rows, {len(df.columns)} columns)")

    if args.ids:
        if "ID" in df.columns:
            df = df[df["ID"].isin(args.ids)]
            print(f"Filtered to {len(df)} rows matching IDs: {args.ids}")
        else:
            print("WARNING: ID column not found — processing all rows")

    all_results = []

    if "oncokb" in modules:
        all_results.extend(run_oncokb(df, args.token))

    if "vep" in modules:
        all_results.extend(run_vep(df))

    # Combined summary
    total_pass     = sum(1 for r in all_results if r["Result"] == "PASS")
    total_mismatch = sum(1 for r in all_results if r["Result"] == "MISMATCH")
    total_error    = sum(1 for r in all_results
                         if r["Result"].startswith("ERROR"))

    print(f"\n{'#'*60}")
    print("COMBINED SUMMARY")
    print(f"{'#'*60}")
    for mod in sorted(modules):
        mod_results = [r for r in all_results if r.get("Module","").lower() == mod]
        p = sum(1 for r in mod_results if r["Result"] == "PASS")
        m = sum(1 for r in mod_results if r["Result"] == "MISMATCH")
        e = sum(1 for r in mod_results if r["Result"].startswith("ERROR"))
        print(f"  {mod.upper():<8}  PASS={p}  MISMATCH={m}  ERROR={e}")
    print(f"  {'TOTAL':<8}  PASS={total_pass}  MISMATCH={total_mismatch}  ERROR={total_error}")
    print(f"{'#'*60}")

    if args.output:
        out_path = Path(args.output)
        pd.DataFrame(all_results).to_csv(out_path, sep="\t", index=False)
        print(f"\nResults written to: {out_path}")

    if total_mismatch > 0 or total_error > 0:
        print("\nRESULT: FAILED")
        sys.exit(1)
    else:
        print("\nRESULT: ALL CHECKS PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
