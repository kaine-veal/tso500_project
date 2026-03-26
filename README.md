# SMART — Somatic Mutation Annotation and Reporting Tool (Dockerised)

## Overview

SMART packages the entire TSO500 variant annotation and filtration pipeline into a single Docker container, eliminating the need for Apptainer/Singularity, separate conda environments, or manual tool installation.

**What's inside the container:**
- GATK 4.6.0.0 (LiftOver)
- Ensembl VEP 114 (with plugin support)
- bcftools, samtools, tabix, bgzip
- Python 3 with all pipeline dependencies (pandas, cyvcf2, etc.)
- OncoKB Annotator
- All SMART pipeline scripts (oncokb2.0.py, vcf2table.py, post_analysis.py)

## Prerequisites

- **Docker** (or Docker Desktop) installed
- **Reference files** already downloaded (chain file, hg38.fa, VEP cache, SpliceAI, REVEL, ClinVar, CIViC, CancerHotspots, LOEUF)
- **OncoKB API token**

## Directory Structure

```
smart/
├── Dockerfile
├── docker-compose.yml
├── entrypoint.sh
├── requirements.txt
├── .dockerignore
├── scripts/   
│   ├── oncokb2.0.py
│   ├── vcf2table.py
│   └── post_analysis.py
└── data/                
    ├── OriginalVcf/
    │   └── *.vcf.gz
    └── TSO500_transcripts_list.txt
```

## Setup

### 1. Build the Docker image

```bash
docker build -t smart:latest .
```

This will take a few minutes on the first build (downloading GATK, VEP, Python packages).

### 4. Edit volume paths

Open `docker-compose.yml` and update the left-hand side of the volume mounts:

```yaml
volumes:
  - ./data:/data                      # Your working directory
  - /mnt/data1/vep_test:/refs:ro      # ← Change to YOUR ref files path
```

## Running SMART

### Option A: Using docker compose

```bash
export ONCOKB_TOKEN=your_token_here

docker compose run --rm smart
```

### Option B: Direct docker run

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v /mnt/data1/vep_test:/refs:ro \
  smart:latest \
  "$ONCOKB_TOKEN" \
  --transcripts-file /data/TSO500_transcripts_list.txt \
  --ref-dir /refs
```

### Common Options

```bash
# Skip PASS filtering (data already filtered)
  --no-pass

# Skip liftover (data already on hg38)
  --no-liftover

# Keep intermediate files for debugging
  --keep-tmp

# Keep per-sample table files
  --keep-tables
```

### Run in background (recommended for large batches)

```bash
docker compose run --rm smart > smart.log 2>&1 &
```

## Output

Results are written inside the mounted `/data` directory:

| Path | Contents |
|------|----------|
| `data/variant_counts.txt` | Per-sample variant counts at each stage |
| `data/Output_Results/` | Final merged results from post_analysis.py |

(If `--keep-tmp` is used, intermediate VCFs are also retained in their respective directories.)

## Key Changes from the Original Pipeline

| Original | SMART (Dockerised) |
|----------|-----------|
| `apptainer exec docker://broadinstitute/gatk:4.6.0.0 /gatk/gatk LiftoverVcf ...` | `gatk LiftoverVcf ...` |
| `apptainer exec docker://ensemblorg/ensembl-vep:release_114.0 vep ...` | `vep ...` |
| `conda env create -f TSO500_environment.yml` | Built into the image |
| `python3 ../oncokb-annotator/MafAnnotator.py` | `python3 /opt/oncokb-annotator/MafAnnotator.py` |

## Troubleshooting

**"Reference resource not found"** — Check that your `--ref-dir` volume mount is correct and the expected subdirectory structure exists.

**Out of memory** — Adjust the memory limit in `docker-compose.yml` under `deploy.resources.limits.memory`.

**VEP cache missing** — The VEP cache (homo_sapiens, 114, GRCh38) must be in your ref-dir. Download it with:
```bash
cd /path/to/your/refs
mkdir -p homo_sapiens
# Download from Ensembl FTP
```

**Permission errors** — Ensure the data directory is writable. You may need to run:
```bash
docker run --rm -u $(id -u):$(id -g) ...
```
