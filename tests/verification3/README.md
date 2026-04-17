# SMART Pipeline — Verification 3: Transcript Prioritisation Impact

## Purpose

This verification demonstrates that the preferred-transcript whitelist meaningfully
changes the VEP annotation output, and that the pipeline consistently selects the
correct transcript for clinically actionable variants.

The pipeline is run twice on the same VCF using two different transcript files.
`verify.py` then confirms:

1. **DIFFER** — targeted variants produce different VEP annotations between runs
   (the transcript whitelist is actually driving annotation selection)
2. **PASS** — each run's annotation is confirmed correct by Ensembl REST API
   using that run's expected transcript

## Test Cases

| Variant | Transcript A (preferred) | Result A | Transcript B (alternative) | Result B |
|---|---|---|---|---|
| EGFR L858R (chr7:55191822 T>G) | NM_005228.5 | p.Leu858Arg (**LEVEL_1** in OncoKB) | NM_001346897.2 | p.Leu813Arg (unknown in OncoKB) |
| TP53 R175H (chr17:7675088 C>T) | NM_000546.6 | p.Arg175His (hotspot, R175H) | NM_001126115.2 | p.Arg43His (different residue entirely) |
| CDKN2A/CDKN2B DEL (chr9:21990000-22005000) | NM_058195.4 | Gene = **CDKN2A** | NM_004936.4 | Gene = **CDKN2B** (different gene reported) |

### Why this matters

- **EGFR L858R**: OncoKB annotation depends entirely on the protein change sent.
  NM_005228 gives the canonical p.Leu858Arg (FDA-approved targeted therapy at LEVEL_1).
  A wrong transcript gives p.Leu813Arg — an unrecognised variant with no therapeutic
  implication, even though the genomic position is identical.

- **TP53 R175H**: Multiple TP53 isoforms exist with different N-terminal lengths.
  NM_001126115 encodes a shorter isoform; the same codon change lands at position 43
  instead of 175, making the hotspot annotation entirely incorrect.

- **CDKN2A/CDKN2B deletion**: The deletion at chr9:21990000-22005000 spans a region
  that overlaps both CDKN2A (p14ARF isoform, NM_058195) and CDKN2B (NM_004936).
  Depending on which transcript is preferred, the pipeline reports a different gene
  as the affected target — a clinically significant distinction.

## Transcript files

| File | Contents |
|---|---|
| `transcripts_A.txt` | NM_005228.5 (EGFR), NM_000546.6 (TP53), NM_058195.4 (CDKN2A p14ARF) |
| `transcripts_B.txt` | NM_001346897.2 (EGFR alt), NM_001126115.2 (TP53 short isoform), NM_004936.4 (CDKN2B) |

## How to run

```bash
export ONCOKB_TOKEN="<your-token>"
bash tests/verification3/run_verification3.sh
```

This will:
1. Run the pipeline with `transcripts_A.txt` → output in `output_A/`
2. Run the pipeline with `transcripts_B.txt` → output in `output_B/`
3. Run `verify.py` comparing both MAF outputs

## verify.py output

```
DIFFER  — annotations differ between runs (expected for targeted variants)
SAME    — annotations are identical between runs (unexpected — investigate)
PASS    — pipeline annotation matches Ensembl REST for the expected transcript
MISMATCH — pipeline annotation does not match Ensembl REST
ERROR   — variant missing from MAF or API error
```

Results are written to `results.tsv`.

Exit code `0` = all checks passed. Exit code `1` = any SAME, MISMATCH, or ERROR.
