# Clinical Variant Interpretation & Reporting Engine (Phase 2)

Extends the [NGS variant calling pipeline](../README.md) with dual-framework clinical variant classification — germline (ACMG/AMP 2015) and somatic (AMP/ASCO/CAP 2017) — plus an automated HTML/PDF clinical report generator.

## Status

Germline arm functional end-to-end:
1. ACMG/AMP 2015 combining-rules classifier — 23/23 tests passing against Richards et al. 2015 Table 5
2. ClinVar evidence extraction (PP5/BP6)
3. gnomAD population frequency evidence (BA1/BS1/PM2) via public API with retry/backoff
4. PVS1 (null variant), gated by a curated gene-mechanism list (JAG1, ASXL1; GNAS deliberately excluded)
5. HTML report generator (Jinja2)

Somatic arm (AMP/ASCO/CAP 2017 tiering) not yet started.

## Data source risk assessment

Verified before use, given repeated large-file failures in Phase 1:

| Source | Size/access | Decision |
|---|---|---|
| gnomAD (full) | Hundreds of GB | Avoided; per-variant public API used instead |
| dbNSFP | Tens of GB | Deferred; PP3/BP4 not yet automated |
| ClinVar VCF | 184MB, plain download | Used directly |
| COSMIC | Free academic, requires registration+SFTP | Deferred to somatic phase |
| OncoKB | Free academic, API-only, no bulk download | Deferred to somatic phase |

## Running

```bash
cd clinical-interpretation
pip install pytest requests jinja2

# Tests
python3 -m pytest tests/ -v

# Evidence data
bash scripts/download_clinvar.sh

# Classify
python3 germline/classify_variants.py \
    --vcf ../results/annotation/NA12878.panel_filtered.annotated.vcf.gz \
    --clinvar resources/clinvar_chr20.vcf.gz \
    --output results/classifications.csv \
    --query-gnomad

# Report
python3 reports/generate_report.py \
    --csv results/classifications.csv \
    --sample-id NA12878 \
    --output reports/NA12878_germline_report.html
```

## References

- Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med.* 2015;17(5):405-424.
- Li MM, Datto M, Duncavage EJ, et al. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer. *J Mol Diagn.* 2017;19(1):4-23.
