# SMART — Somatic Mutation Annotation and Reporting Tool

<p align="center"><em>A Dockerised pipeline for TSO500 variant annotation, filtration, and clinical reporting.</em></p>

---

## Overview

SMART automates the end-to-end processing of VCF files from the Illumina TSO500 panel. Starting from raw Dragen VCFs on hg19/GRCh37, it performs PASS filtering, coordinate liftover to hg38, comprehensive functional annotation via Ensembl VEP, clinical annotation via OncoKB, and produces a final analysis-ready table with transcript-prioritised results.

Everything runs inside a single Docker container — no Apptainer, Singularity, or conda environment required.

---

## Workflow

```
 ┌──────────────┐
 │  Input VCFs  │   Dragen .vcf.gz files (hg19)
 │  OriginalVcf │
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 1. PASS      │   Keep only PASS-filtered variants
 │    Filter    │   (awk + bgzip + tabix)
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 2. LiftOver  │   hg19 → hg38 coordinate conversion
 │    (GATK)    │   Rejected variants saved separately
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 3. VEP       │   Functional annotation with plugins:
 │  Annotation  │   SpliceAI, REVEL, ClinVar, CIViC,
 │              │   CancerHotspots, LOEUF
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 4. OncoKB    │   Clinical annotation via API:
 │  Annotation  │   Oncogenicity, treatment levels,
 │              │   diagnostic & prognostic implications
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 5. VCF →     │   Transcript-prioritised table with
 │    Table     │   all VEP + OncoKB fields
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 6. MAF       │   OncoKB MafAnnotator for
 │  Annotation  │   standardised MAF output
 └──────┬───────┘
        │
        ▼
 ┌──────────────┐
 │ 7. Post      │   Merge all sample tables into
 │  Analysis    │   final combined results
 └──────────────┘
```

---

## Transcript Selection Strategy

A key design principle of SMART is consistent transcript prioritisation. Both the OncoKB annotator and the table generator use the same 3-tier logic to ensure the transcript shown in the final output matches the one used to query OncoKB:

**Tier 1 — Preferred List:** Does the annotation contain a RefSeq NM ID (version-agnostic) found in the TSO500 transcript whitelist?

**Tier 2 — MANE Select / MANE Plus Clinical:** If no Tier 1 match, use the transcript tagged as MANE Select or MANE Plus Clinical by Ensembl/NCBI.

**Tier 3 — Fallback:** If neither applies, use the first transcript reported by VEP.

---

## What's Inside the Container

| Tool | Version | Purpose |
|------|---------|---------|
| GATK | 4.6.0.0 | LiftOver (hg19 → hg38) |
| Ensembl VEP | 114.0 | Functional annotation |
| SpliceAI plugin | 1.3 | Splice-site impact prediction |
| REVEL plugin | 1.3 | Missense pathogenicity scoring |
| bcftools / samtools / tabix | System | VCF manipulation |
| Python 3.10 | + pandas, cyvcf2, requests | Pipeline scripts |
| OncoKB Annotator | Latest | MafAnnotator for MAF output |

---

## Prerequisites

- **Docker** (v20.10+ or Docker Desktop)
- **OncoKB API token** — obtain from [oncokb.org](https://www.oncokb.org)
- **Reference files** downloaded and organised (see below)

---

## Reference Files Setup

SMART expects a reference directory with the following structure:

```
/path/to/refs/
├── liftover/
│   ├── hg19ToHg38.over.chain       # UCSC chain file
│   └── hg38.fa                      # GRCh38 reference genome + .fai + .dict
├── Plugins/                         # VEP plugin files (SpliceAI.pm, REVEL.pm, LOEUF.pm)
├── SpliceAI/
│   ├── spliceai_scores.raw.snv.hg38.vcf.gz      (+.tbi)
│   └── spliceai_scores.raw.indel.hg38.vcf.gz    (+.tbi)
├── REVEL/
│   └── new_tabbed_revel_grch38.tsv.gz            (+.tbi)
├── ClinVar/
│   └── clinvar.vcf.gz                            (+.tbi)
├── CIVIC/
│   └── civic_01_10_25.vcf.gz                     (+.tbi)
├── gnomAD_constraints/
│   └── loeuf_dataset_grch38.tsv.gz               (+.tbi)
├── CancerHotSpots/
│   └── hg38.hotspots_changv2_gao_nc.vcf.gz      (+.tbi)
└── homo_sapiens/                    # VEP cache directory
    └── 114_GRCh38/
```

**Download sources:**

| Resource | Source |
|----------|--------|
| Chain file | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` |
| hg38 reference | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` |
| VEP cache | Ensembl FTP (`homo_sapiens_vep_114_GRCh38.tar.gz`) |
| SpliceAI | Illumina basespace (requires account) |
| ClinVar | NCBI FTP |
| CIViC | CIViC downloads page |

---

## Quick Start

### 1. Clone and prepare

```bash
git clone https://github.com/kaine-veal/tso500_project.git
cd tso500_project

# Create the Docker build context
mkdir -p docker/scripts
cp oncokb2.0.py vcf2table.py post_analysis.py docker/scripts/
cp Dockerfile docker-compose.yml entrypoint.sh requirements.txt .dockerignore docker/
```

### 2. Add your data

```bash
mkdir -p docker/data/OriginalVcf
cp /path/to/your/*.vcf.gz docker/data/OriginalVcf/
cp TSO500_transcripts_list.txt docker/data/
```

### 3. Build

```bash
cd docker
docker build -t smart:latest .
```

### 4. Run

```bash
export ONCOKB_TOKEN=your_token_here

docker run --rm \
  -v $(pwd)/data:/data \
  -v /path/to/your/refs:/refs:ro \
  smart:latest \
  "$ONCOKB_TOKEN" \
  --transcripts-file /data/TSO500_transcripts_list.txt \
  --ref-dir /refs
```

Or with docker compose (edit `docker-compose.yml` volume paths first):

```bash
docker compose run --rm smart
```

### 5. Run in background (recommended for large batches)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v /path/to/your/refs:/refs:ro \
  smart:latest \
  "$ONCOKB_TOKEN" \
  --transcripts-file /data/TSO500_transcripts_list.txt \
  --ref-dir /refs \
  > smart.log 2>&1 &
```

---

## Command-Line Options

```
Usage:
  smart <ONCOKB_TOKEN> --transcripts-file FILE --ref-dir DIR [OPTIONS]

Required:
  <ONCOKB_TOKEN>                   OncoKB API token
  --transcripts-file FILE          Transcript whitelist for prioritisation
  --ref-dir DIR                    Reference resources directory

Options:
  --pass / --no-pass               Enable/disable PASS filtering (default: ON)
  --liftover / --no-liftover       Enable/disable hg19→hg38 liftover (default: ON)
  --clean-tmp / --keep-tmp         Delete/keep intermediate files (default: clean)
  --clean-tables / --keep-tables   Delete/keep per-sample tables after
                                   post_analysis merging (default: clean)
  --help                           Show help
```

---

## Output

All results are written to the mounted `/data` directory:

| File | Description |
|------|-------------|
| `variant_counts.txt` | Per-sample variant counts at each pipeline stage |
| `Output_Results/` | Final merged results from post-analysis |

When `--keep-tmp` is used, intermediate files are retained:

| Directory | Contents |
|-----------|----------|
| `FilteredVcf/` | PASS-only VCFs |
| `LiftOverVcf/` | hg38-lifted VCFs |
| `LiftOverVcf/Rejected/` | Variants that failed liftover |
| `AnnotatedVcf/` | VEP-annotated VCFs |
| `OncoKB_VCF/` | OncoKB-annotated VCFs |
| `Table/` | Per-sample CSV tables |
| `FINAL_Table/` | Per-sample MAF files |

---

## Output Annotations

The final table includes annotations from multiple sources. Key field groups:

**Core VCF:** CHROM, POS, ID, REF, ALT, QUAL, FILTER, FORMAT

**VEP Functional:** Consequence, IMPACT, SYMBOL, HGVSc, HGVSp, BIOTYPE, EXON, INTRON, SIFT, PolyPhen, VARIANT_CLASS, CANONICAL, MANE_SELECT

**Population Frequency:** gnomADe/gnomADg allele frequencies across populations, MAX_AF

**Splice & Pathogenicity:** SpliceAI delta scores (AG/AL/DG/DL), REVEL missense score, LOEUF constraint

**Clinical Databases:** ClinVar (significance, review status, conditions), CIViC (variant type, consequence), CancerHotspots (hotspot and 3D hotspot flags)

**OncoKB Clinical:** Oncogenic classification, mutation effect, therapeutic/diagnostic/prognostic levels, treatment summaries, gene/variant summaries

---

## Variant Classification

SMART classifies variants to route them to the correct OncoKB API endpoint:

| Variant Type | Examples | OncoKB Endpoint |
|-------------|----------|-----------------|
| Mutations | SNVs, indels, MantaINS, MantaBND | `/annotate/mutations/byProteinChange` |
| Copy Number | MantaDUP, MantaDEL, GAIN, LOSS | `/annotate/copyNumberAlterations` |

Structural insertions (MantaINS) are intentionally routed to the mutation endpoint because they alter gene sequence and produce specific protein changes (e.g. EGFR Exon 20 insertions) that require distinct targeted therapies, unlike simple gene dosage changes.

---

## Tumor Type Handling

SMART supports two modes for OncoKB tumour type context:

**`--tumor_mode filename`** (default in standalone script): Parses the input filename against a built-in cancer type dictionary to infer the tumour type for context-specific evidence levels.

**`--tumor_mode generic`** (default in Docker entrypoint): Uses "UNKNOWN" for pan-cancer annotation. Recommended per OncoKB documentation for broader evidence capture.

Supported filename-inferred types: brain, breast, cholangiocarcinoma, colon, endometrial, HNSC, SNUC, lung, melanoma, ovarian, pancreatic, prostate, renal, sarcoma, thyroid.

---

## Tool Versions

| Component | Version |
|-----------|---------|
| VCF format | 4.2 |
| GATK LiftoverVcf | 4.6.0.0 |
| Ensembl VEP | 114.2 |
| OncoKB | 5.4 |
| SpliceAI | 1.3 |
| REVEL | 1.3 |
| CIViC | 01_10_25 |
| Reference genome | GRCh38 (Homo sapiens 114) |

---

## Troubleshooting

**"Required resource not found"**
Check that your `--ref-dir` volume mount is correct and the expected subdirectory structure exists. All `.vcf.gz` reference files need accompanying `.tbi` index files.

**Out of memory**
VEP annotation can be memory-intensive. Adjust the memory limit in `docker-compose.yml` under `deploy.resources.limits.memory`, or pass `--memory=32g` to `docker run`.

**VEP cache missing**
The VEP cache (`homo_sapiens/114_GRCh38/`) must exist inside your ref-dir. Download from Ensembl FTP and extract.

**Permission errors on output files**
The container runs as root by default. To match your host user:
```bash
docker run --rm -u $(id -u):$(id -g) ...
```

**OncoKB API errors**
Verify your token is valid and has not expired. The pipeline logs API errors to stderr — check your log file for `OncoKB mutation API ERROR` messages.

**Empty VCF after PASS filtering**
Some samples may have zero PASS variants. The pipeline will process these without error but the sample will show 0 counts in `variant_counts.txt`.

---

## Authors

Manuel, Mani, Kaine, and Ian

---

## License

See repository for licence details.
