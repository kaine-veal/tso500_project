# SMART — Somatics Annotation Reported Tool

**TSO500 variant annotation and reporting pipeline**

[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://www.python.org/)
[![Shell](https://img.shields.io/badge/Shell-Bash-green)](https://www.gnu.org/software/bash/)
[![VEP](https://img.shields.io/badge/VEP-114.2-orange)](https://www.ensembl.org/vep)
[![OncoKB](https://img.shields.io/badge/OncoKB-v5.4-red)](https://www.oncokb.org/)

SMART is an end-to-end somatic variant annotation pipeline for **Illumina TSO500** panel data. Starting from Dragen-generated VCF files aligned to hg19/GRCh37, it lifts coordinates to GRCh38, applies comprehensive functional annotation via VEP, enriches variants with clinical evidence from OncoKB, CIViC, ClinVar and Cancer Hotspots, and delivers three audience-targeted output tables ready for bioinformatic analysis and clinical review.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Step-by-Step Details](#step-by-step-details)
- [Output Files](#output-files)
- [Annotation Fields Reference](#annotation-fields-reference)
- [Tool Versions](#tool-versions)
- [Data Location](#data-location)
- [Known Issues and TODO](#known-issues-and-todo)
- [Acknowledgements](#acknowledgements)

---

## Overview

The pipeline processes Illumina TSO500 VCF files through six sequential stages:

1. **PASS filtering** — retains only high-confidence variant calls
2. **LiftOver** — converts coordinates from hg19 to hg38 using GATK
3. **VEP annotation** — adds functional, population, splicing, and clinical annotations
4. **OncoKB annotation** — enriches variants with therapeutic, diagnostic and prognostic evidence
5. **Table generation** — converts annotated VCFs to structured tabular output with consistent transcript selection
6. **Post-analysis** — merges all sample tables and produces the final tiered output files

---

## Pipeline Architecture

```
OriginalVcf/
  └── *.vcf.gz  (Dragen TSO500 output, hg19)
        │
        ▼  Step 1: PASS filter (awk + bgzip)
  FilteredVcf/
        │
        ▼  Step 2: LiftOver hg19 → hg38 (GATK LiftoverVcf via Apptainer)
  LiftOverVcf/
    └── Rejected/   (variants that failed coordinate conversion)
        │
        ▼  Step 3: VEP annotation (VEP 114.2 + plugins)
  AnnotatedVcf/
        │
        ▼  Step 4: OncoKB annotation (oncokb2.0.py)
  OncoKB_VCF/
        │
        ▼  Step 5: VCF → table (vcf2table.py)
  FINAL_Table/
        │
        ▼  Step 6: Post-analysis (post_analysis.py + build_config.py)
  Output_Results/
    ├── Final_result_tier1.maf     — all non-redundant fields, single header (MAF format)
    ├── Final_result_tier2.tsv     — ~400 fields for bioinformaticians, two-row header
    └── Final_result_tier3.tsv     — ~200 fields for clinical scientists, two-row header
```

---

## Prerequisites

### System dependencies

| Tool | Purpose | Notes |
|---|---|---|
| Apptainer / Singularity | Containerised execution of GATK | Requires system admin installation |
| GATK ≥ 4.6.0.0 | LiftOver | Run via Apptainer image |
| Ensembl VEP 114.2 | Functional annotation | [VEP install docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) |
| bcftools | VCF manipulation | Via Apptainer or conda |
| Conda / Mamba | Python environment management | [Miniforge](https://github.com/conda-forge/miniforge) |

### VEP plugins and custom databases

The following must be installed and configured in your VEP data directory before running:

| Plugin / Database | Version | Purpose |
|---|---|---|
| SpliceAI | 1.3 | Splice impact prediction |
| REVEL | v1.3 | Missense pathogenicity scores |
| ClinVar | ClinVar_2024-12 or later | Germline and somatic clinical classifications |
| CIViC | civic_01_10_25 or later | Clinical variant interpretations |
| Cancer Hotspots | changv2 | Recurrent somatic mutation hotspots |

### Reference files

| File | Source | Download URL |
|---|---|---|
| `hg19ToHg38.over.chain` | UCSC | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` |
| `hg38.fa` | NCBI GRCh38 analysis set | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` |

### OncoKB API token

A valid OncoKB API token is required. Register at [oncokb.org](https://www.oncokb.org/) to obtain one.

---

## Repository Structure

```
tso500_project/
├── DoIt.sh                      # Main pipeline entry point
├── main.sh                      # Core pipeline logic (called by DoIt.sh)
├── get_ref_files.sh             # Helper to download reference files
├── oncokb2.0.py                 # OncoKB annotation engine
├── vcf2table.py                 # VCF → tabular conversion with transcript selection
├── post_analysis.py             # Merge, tier, and produce final output files
├── ConfigORIGINAL.yaml          # Human-maintained field metadata template
├── create_config/
│   └── build_config.py          # Generates Config.yaml from ConfigORIGINAL.yaml
├── TSO500_environment.yml       # Conda environment specification
├── TSO500_transcripts_list.txt  # Preferred transcript list (NM_ IDs, no version numbers)
├── Files4ThisProject/           # Supporting reference files
├── OriginalVcf/                 # Input: Dragen TSO500 VCF files (*.vcf.gz)
└── Output_Results/              # Final output tables
```

---

## Quick Start

### 1. Create the Conda environment

```bash
conda env create -f TSO500_environment.yml -n TSO500
conda activate TSO500
```

### 2. Place input files

Put your Dragen VCF files into `./OriginalVcf/`. Files must follow this naming convention:

```
D21_11111_22222_TumourType_33333_S1.dragen.concat_snv_sv.vcf.gz
```

The tumour type field in the filename is parsed by `oncokb2.0.py` to select cancer-specific OncoKB evidence.

### 3. Configure paths

Edit the top of `main.sh` to point to:
- `CHAIN` — path to `hg19ToHg38.over.chain`
- `REF` — path to `hg38.fa`
- VEP installation directory, cache, plugins, and custom database paths
- Apptainer image for GATK

### 4. Run the pipeline

```bash
./DoIt.sh YOUR_ONCOKB_TOKEN
```

For large cohorts, run in the background:

```bash
nohup ./DoIt.sh YOUR_ONCOKB_TOKEN > DoIt.log 2>&1 &
```

Progress is written to stdout/log in real time. Processing ~420 samples typically takes several hours depending on server resources.

---

## Step-by-Step Details

### Step 1: PASS Filtering

Retains only variants where the `FILTER` column is exactly `PASS`, using `awk` with `zcat` and `bgzip`. All subsequent steps work exclusively on high-confidence calls.

### Step 2: LiftOver (hg19 → hg38)

Runs **GATK `LiftoverVcf`** via Apptainer to convert genomic coordinates from hg19 (GRCh37) to hg38 (GRCh38). Variants that cannot be converted are written to `LiftOverVcf/Rejected/` for inspection.

### Step 3: VEP Annotation

Runs **Ensembl VEP 114.2** with `--everything` plus the following custom annotations:

| Plugin / Database | Adds |
|---|---|
| SpliceAI (v1.3) | Splice acceptor/donor gain and loss delta scores and positions |
| REVEL (v1.3) | Missense pathogenicity composite score |
| ClinVar | Germline, oncogenic, and somatic clinical impact classifications |
| CIViC | Clinical variant interpretations and evidence levels |
| Cancer Hotspots | Recurrent somatic hotspot flags (protein-level, 3D, and non-coding) |
| gnomAD exomes + genomes | Per-population allele frequencies |
| LOEUF | Gene-level loss-of-function intolerance |

### Step 4: OncoKB Annotation (`oncokb2.0.py`)

The core clinical annotation engine. Queries the **OncoKB API** to retrieve therapeutic, diagnostic, and prognostic evidence for each variant.

#### Tumour type inference

The script parses the VCF filename against a predefined `CANCER_MAP` dictionary to detect the tumour type and pass it to OncoKB for cancer-specific evidence levels (e.g. Level 1 for EGFR L858R in Lung Cancer vs. Level 3 in other cancer types). Use `--tumor_mode=generic` to force pan-cancer queries instead.

#### Variant routing

| Variant type | OncoKB endpoint |
|---|---|
| SNV, Indel, `MantaINS` | `/annotate/mutations/byProteinChange` |
| `MantaDUP`, `MantaDEL` | `/annotate/copyNumberAlterations` |

> **Note on `MantaINS` routing:** Structural insertions are intentionally sent to the mutation endpoint rather than the CNV endpoint. Insertions such as EGFR Exon 20 insertions alter internal gene sequence and require specific targeted therapies entirely distinct from gene amplification treatments. The OncoKB CNV endpoint also does not accept `INSERTION` as a valid alteration type.

#### Transcript selection — 3-tier priority system

To ensure the clinically relevant isoform is used for OncoKB queries, the script applies the following priority hierarchy to the VEP `CSQ` field:

| Tier | Rule |
|---|---|
| **Tier 1 (Preferred)** | RefSeq NM_ ID matching `TSO500_transcripts_list.txt` (version-agnostic) |
| **Tier 2 (MANE)** | Transcript tagged as MANE Select or MANE Plus Clinical |
| **Tier 3 (Fallback)** | First transcript reported by VEP |

### Step 5: VCF to Table (`vcf2table.py`)

Converts the fully annotated VCF to a structured tabular file. Applies the **identical 3-tier transcript selection logic** as Step 4 to guarantee that the transcript shown in the output table is the same one used to query OncoKB.

### Step 6: Post-Analysis (`post_analysis.py`)

Merges all per-sample tables and produces three audience-targeted output files. Field metadata (descriptions, sources, versions, and tier assignments) is read from `Config.yaml`, which is generated by running:

```bash
python create_config/build_config.py
```

This script reads `ConfigORIGINAL.yaml` (the human-maintained field metadata template) and outputs `Config.yaml`, which is consumed by `post_analysis.py` at runtime. Edit `ConfigORIGINAL.yaml` to update field descriptions, sources, versions, or tier assignments — never edit `Config.yaml` directly.

---

## Output Files

Three files are produced in `Output_Results/`:

| File | Audience | Approx. columns | Header |
|---|---|---|---|
| `Final_result_tier1.maf` | All users / downstream tools | ~900+ (all non-redundant) | Single row (field names only) |
| `Final_result_tier2.tsv` | Bioinformaticians | ~400 | Two rows: field names + `description \| source \| version` |
| `Final_result_tier3.tsv` | Clinical scientists | ~200 | Two rows: field names + `description \| source \| version` |

> ⚠️ `Final_result_tier1.maf` uses the `.maf` extension but **has not yet been validated** for compatibility with MAF-consuming tools (e.g. cBioPortal, vcf2maf). Verify that all required mandatory MAF columns are present before using with downstream tools.

### Tier definitions

Every field appears in exactly one tier. Higher-numbered tiers are strict subsets of lower-numbered ones — any field in Tier 3 is also in Tier 2 and Tier 1.

| Tier | Appears in |
|---|---|
| **T3** | All three outputs (MAF + bioinformatician TSV + clinician TSV) |
| **T2** | MAF + bioinformatician TSV only |
| **T1** | MAF only |
| **DROP** | Excluded from all outputs |

### Why fields are dropped

| Reason | Examples |
|---|---|
| Superseded by expanded columns | `ONCOKB_JSON`, `ONCOKB_treatments`, `ONCOKB_diagnosticImplications` → content fully expanded into `ONCOKB_TX_*` and `ONCOKB_DIAG_*` |
| Internal query metadata | `ONCOKB_query.referenceGenome`, `ONCOKB_query.hugoSymbol` — pipeline input parameters |
| Truncated legacy names | `ONCOKB_genefx`, `ONCOKB_oncog`, `ONCOKB_highes` — superseded by properly named columns |
| Raw unexpanded annotation | `CIViC_CSQ` — fully parsed into `CIViC_CSQ_*` subfields |
| Duplicate column | `CHROM` — exact duplicate of the MAF-standard `Chromosome` column |

---

## Annotation Fields Reference

### Core VCF

| Field | Description |
|---|---|
| `Chromosome` | Chromosome name |
| `Start_Position` | Variant start position |
| `End_Position` | Variant end position |
| `REF` | Reference allele |
| `ALT` | Alternate allele |
| `FILTER` | Filter status (all `PASS` after Step 1) |
| `Existing_variation` | Known variant identifiers (rsID, COSV) |

### Genotype / FORMAT

| Field | Description |
|---|---|
| `GT` | Genotype call |
| `AD` | Allelic depths for ref and alt alleles |
| `VAF` | Variant allele fraction |
| `DP` | Read depth |
| `F1R2` / `F2R1` | Reads supporting allele by orientation |
| `SB` | Strand bias (Fisher's Exact Test components) |
| `SQ` | Somatic quality score |

### VEP

| Field | Description |
|---|---|
| `SYMBOL` | Gene symbol |
| `NM_Transcript` | Selected RefSeq transcript |
| `HGVSp_Short` | Abbreviated protein change (e.g. p.V600E) |
| `HGVSc` | HGVS coding DNA notation |
| `Consequence` | Sequence Ontology consequence |
| `IMPACT` | HIGH / MODERATE / LOW / MODIFIER |
| `EXON` / `INTRON` | Exon or intron number |
| `MANE_SELECT` | MANE Select transcript |
| `SIFT` | SIFT pathogenicity prediction |
| `PolyPhen` | PolyPhen pathogenicity prediction |
| `REVEL` | Missense pathogenicity composite score |
| `LOEUF` | Loss-of-function intolerance score |

### Population Frequency

| Field | Description |
|---|---|
| `AF` | Global allele frequency (1000 Genomes) |
| `gnomADe_AF` | gnomAD exomes global AF |
| `gnomADg_AF` | gnomAD genomes global AF |
| `MAX_AF` / `MAX_AF_POPS` | Maximum AF and population where it was observed |
| `gnomADe_AFR_AF` … `gnomADe_SAS_AF` | Per-population gnomAD exomes AFs |
| `gnomADg_AFR_AF` … `gnomADg_SAS_AF` | Per-population gnomAD genomes AFs |

### SpliceAI (v1.3)

| Field | Description |
|---|---|
| `SpliceAI_pred_DS_AG` | Delta score: acceptor gain (0–1) |
| `SpliceAI_pred_DS_AL` | Delta score: acceptor loss (0–1) |
| `SpliceAI_pred_DS_DG` | Delta score: donor gain (0–1) |
| `SpliceAI_pred_DS_DL` | Delta score: donor loss (0–1) |
| `SpliceAI_pred_DP_*` | Distance to predicted affected splice site |
| `SpliceAI_cutoff` | PASS/FAIL flag relative to pipeline threshold |

### ClinVar

| Field | Description |
|---|---|
| `ClinVar_CLNSIG` | Aggregate germline classification |
| `ClinVar_CLNREVSTAT` | Germline review status |
| `ClinVar_CLNDN` | Associated disease name |
| `ClinVar_ONC` | Oncogenicity classification |
| `ClinVar_ONCREVSTAT` | Oncogenicity review status |
| `ClinVar_SCI` | Somatic clinical impact |
| `ClinVar_SCIREVSTAT` | Somatic clinical impact review status |
| `ClinVar_ORIGIN` | Allele origin (0=unknown; 1=germline; 2=somatic; etc.) |

### Cancer Hotspots

| Field | Description |
|---|---|
| `CancerHotspots` | Hotspot overlap flag |
| `CancerHotspots_HOTSPOT` | Protein-level hotspot annotation |
| `CancerHotspots_HOTSPOT3D` | 3D structural cluster hotspot flag |
| `CancerHotspots_HOTSPOTNC` | Non-coding hotspot flag |

### CIViC

| Field | Description |
|---|---|
| `CIViC_GN` | Gene symbol |
| `CIViC_VT` | Variant name |
| `CIViC_CSQ_CIViC Evidence Level` | Evidence level (A–E) |
| `CIViC_CSQ_CIViC Entity Significance` | Clinical significance |
| `CIViC_CSQ_CIViC Entity Therapies` | Associated therapies |
| `CIViC_CSQ_CIViC Assertion AMP Category` | AMP/ASCO/CAP classification |

### OncoKB

| Field | Description |
|---|---|
| `ANNOTATED` | Successfully annotated by OncoKB (TRUE/FALSE) |
| `GENE_IN_ONCOKB` | Gene curated in OncoKB |
| `VARIANT_IN_ONCOKB` | Variant curated in OncoKB |
| `ONCOGENIC` | Oncogenic / Likely Oncogenic / Likely Neutral / Unknown |
| `MUTATION_EFFECT` | Gain-of-function / Loss-of-function / Neutral / etc. |
| `ONCOKB_HOTSPOT` | OncoKB hotspot flag |
| `ONCOKB_SENS_LVL` | Highest therapeutic sensitivity level |
| `ONCOKB_FDA_LVL` | Highest FDA evidence level |
| `ONCOKB_DIAG_LVL` | Highest diagnostic implication level |
| `HIGHEST_LEVEL` | Highest overall therapeutic level (R1 > 1 > 2 > 3A > 3B > 4 > R2) |
| `LEVEL_1` … `LEVEL_R2` | Therapies at each evidence level |
| `ONCOKB_TX_N_level` | Evidence level for treatment entry N (N = 0–29) |
| `ONCOKB_TX_N_drugs` | Drug(s) for treatment entry N |
| `ONCOKB_TX_N_description` | Evidence narrative for treatment entry N |
| `ONCOKB_DIAG_N_levelOfEvidence` | Diagnostic level for entry N (N = 0–3) |
| `ONCOKB_DIAG_N_tumorType.name` | Tumour type for diagnostic entry N |
| `ONCOKB_DATA_VERSION` | OncoKB dataset version |
| `ONCOKB_LAST_UPDATE` | Date of last OncoKB record update |

---

## Tool Versions

| Tool | Version |
|---|---|
| VCF format | 4.2 |
| GATK LiftoverVcf | 4.6.0.0 |
| Ensembl VEP | 114.2 |
| Homo sapiens cache | 114 GRCh38 |
| OncoKB annotator | 5.4 |
| SpliceAI | 1.3 |
| REVEL | v1.3 |
| CIViC | civic_01_10_25 |
| ClinVar | ClinVar_2024-12 |
| Cancer Hotspots | changv2 |
| Python | 3.10 |

---

## Data Location

Processed data on the archive server:

```
/mnt/ArchiveData/results/TSO500_WESSEX/
├── excel_and_other_files_provided/   # Documents shared by team / Birmingham Bioinformatics
├── raw/
│   └── files/                        # All TSO500 FASTQ (.ora) files
├── README.md                         # Data organisation documentation
├── scripts/                          # Integrity checks and sample grouping scripts
└── VCFs/
    └── files/
        ├── *.vcf.gz                  # All TSO500 VCF files
        └── AnnotatedVEP/             # 418 PASS, lifted, VEP-annotated VCFs + DoIt.log
```

---

## Known Issues and TODO

- [ ] End-to-end reanalysis of all 420 samples to verify no variants are missed at any stage
- [ ] Investigate LiftOver rejected variants — expected to be very few since only PASS variants reach that step
- [ ] Compare tumour-specific vs. generic OncoKB query modes and determine which provides better clinical annotation coverage
- [ ] Validation: create a synthetic VCF with known OncoKB-annotated variants (e.g. BRAF V600E, EGFR L858R) to verify the full annotation chain end-to-end
- [ ] Validate `Final_result_tier1.maf` compatibility with MAF-consuming tools (cBioPortal, vcf2maf, GATK) — mandatory MAF columns and format need verification before production use
- [ ] Run VEP in tabular mode as a cross-check on all non-OncoKB annotation fields
- [ ] Confirm team decision on Tier 3 (clinical scientist) field selection
- [ ] Add sample name propagation to the final output tables
- [ ] Update descriptions for `ONCOKB_PROGNOSTIC_IMPLICATIONS`, `ONCOKB_variantSummary`, and `ONCOKB_VUS.1` in `ConfigORIGINAL.yaml`

---

## Acknowledgements

Pipeline developed as part of the TSO500 somatic variant annotation project at University Hospital Southampton NHS Foundation Trust. Clinical annotation powered by [OncoKB](https://www.oncokb.org/), [CIViC](https://civicdb.org/), [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [Cancer Hotspots](https://www.cancerhotspots.org/), [SpliceAI](https://github.com/Illumina/SpliceAI), [REVEL](https://sites.google.com/site/revelgenomics/), and [Ensembl VEP](https://www.ensembl.org/vep).
