# SMART Pipeline — Annotation Verification 1

## Purpose

This verification set checks that the SMART annotation pipeline produces
correct and clinically meaningful results for well-characterised oncogenic
hotspot variants. All variants are in GRCh38 coordinates and have established
annotations in ClinVar, CIViC, and OncoKB.

## Approach

A synthetic VCF (`verification1.vcf`) was created containing 10 known
oncogenic hotspot SNVs and 4 Manta-format CNVs (structural variants).
Each variant was selected because:

1. Its pathogenicity and clinical significance are well established
2. It should be annotated by multiple sources in the pipeline (ClinVar, CIViC,
   OncoKB, CancerHotSpots)
3. Its expected gene, consequence, and oncogenicity are known independently of
   the pipeline

A dedicated transcript list (`verification1_transcripts.txt`) is used instead
of the full TSO500 list to ensure VEP picks the correct canonical transcript
for each gene and that protein changes match expectations.

After running SMART on this VCF, the output is inspected manually to confirm
that each variant receives the expected annotation.

## Variants

| Variant | Gene | Transcript | HGVS protein | Expected significance |
|---|---|---|---|---|
| BRAF_V600E | BRAF | NM_004333.6 | p.Val600Glu | Oncogenic. Targeted therapy in melanoma, CRC, NSCLC |
| KRAS_G12D | KRAS | NM_004985.5 | p.Gly12Asp | Oncogenic. Pan-cancer driver |
| KRAS_G13D | KRAS | NM_004985.5 | p.Gly13Asp | Oncogenic. Pan-cancer driver |
| PIK3CA_H1047R | PIK3CA | NM_006218.4 | p.His1047Arg | Oncogenic. Breast, CRC |
| TP53_R175H | TP53 | NM_000546.6 | p.Arg175His | Oncogenic (gain-of-function). Pan-cancer. Note: minus-strand gene — VCF allele is C>T (genomic+) |
| EGFR_L858R | EGFR | NM_005228.5 | p.Leu858Arg | Oncogenic. NSCLC targeted therapy |
| IDH1_R132H | IDH1 | NM_005896.4 | p.Arg132His | Oncogenic. Glioma, AML. Note: minus-strand gene — VCF allele is C>T (genomic+) |
| NRAS_Q61R | NRAS | NM_002524.5 | p.Gln61Arg | Oncogenic. Melanoma, AML |
| PTEN_R130stop | PTEN | NM_000314.8 | p.Arg130* | Likely oncogenic (tumour suppressor loss). Pan-cancer. Note: minus-strand gene — VCF allele is C>T at chr10:87933147 |
| ERBB2_S310F | ERBB2 | NM_004448.4 | p.Ser310Phe | Oncogenic. Breast, bladder. Note: plus-strand gene — VCF allele is C>T at chr17:39711955 |

## Manta CNV variants

| Variant ID | Gene | SVTYPE | GRCh38 coordinates | Transcript | Expected OncoKB |
|---|---|---|---|---|---|
| MantaDUP:ERBB2_AMP | ERBB2 | DUP | chr17:39,687,914-39,730,426 | NM_004448.4 | Oncogenic, Level 1 (breast cancer — trastuzumab) |
| MantaDUP:MET_AMP | MET | DUP | chr7:116,672,196-116,798,037 | NM_000245.4 | Oncogenic, Level 1 (NSCLC — capmatinib/tepotinib) |
| MantaDUP:CDK4_AMP | CDK4 | DUP | chr12:57,748,512-57,793,240 | NM_000075.4 | Oncogenic, Level 1 (liposarcoma — palbociclib) |
| MantaDEL:CDKN2A_DEL | CDKN2A | DEL | chr9:21,967,753-22,009,812 | NM_000077.5 | Oncogenic (pan-cancer TSG) |

These use real Manta VCF format (`<DUP>`/`<DEL>` ALT, `SVTYPE`/`END`/`SVLEN` INFO fields,
`PR`/`SR` FORMAT fields). The pipeline classifies them as CNVs via `MantaDUP:`/`MantaDEL:` ID
prefix and routes them to the OncoKB `/annotate/copyNumberAlterations` endpoint.

## What to check in the output

For each variant verify:

- **Gene and transcript** — correct gene name and the exact transcript from `verification1_transcripts.txt`
- **Protein change** — matches the expected HGVS protein in the table above
- **ClinVar** — CLNSIG should reflect pathogenic/likely pathogenic status
- **CIViC** — evidence entries present for actionable variants (BRAF, EGFR, KRAS, IDH1)
- **OncoKB** — oncogenicity and therapeutic actionability level assigned
- **CancerHotSpots** — HOTSPOT flag present for known hotspot positions
- **REVEL** — score present for missense variants (not applicable to stop variants)
- **SpliceAI** — score present (expected to be low for coding missense variants)

## How to run

```bash
docker run --rm \
  -v /Users/manolodominguez/tso500_project/tests/verification1:/data \
  quay.io/biocontainers/htslib:1.21--h566b1c6_1 \
  sh -c "bgzip /data/verification1.vcf && tabix -p vcf /data/verification1.vcf.gz"
```

Run SMART:
```bash
docker run --rm \
  -v /Users/manolodominguez/tso500_project/tests/verification1:/data \
  -v /Users/manolodominguez/tso500_project/tests/verification1/output:/output \
  -v /Volumes/ExternalSSD/refs:/refs:ro \
  monkiky/smart:1.0.0 \
  "$ONCOKB_TOKEN" \
  --transcripts-file /data/verification1_transcripts.txt \
  --ref-dir /refs \
  --config /data/Config.yaml \
  --input-dir /data \
  --output-dir /output \
  --no-liftover \
  --keep-tmp 2>&1 | tee /Users/manolodominguez/tso500_project/tests/verification1/output/smart.log
```
