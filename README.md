# Clinical Genomics Pipeline — Variant Calling + ACMG/AMP Interpretation

A two-phase clinical genomics pipeline: **Phase 1** is a reproducible, benchmarked WGS variant calling pipeline built with **Snakemake** and **Docker**, **Phase 2** extends it with an **ACMG/AMP 2015** germline variant classifier and automated clinical report generator.

Two variant callers - **GATK HaplotypeCaller** and **Google DeepVariant** are run in parallel and benchmarked head-to-head against the **Genome in a Bottle (GIAB) NA12878** truth set. The best caller's output is annotated with **VEP**, filtered to a clinical gene panel, cross-referenced against **ClinVar** and **gnomAD** and classified into a defensible ACMG/AMP pathogenicity call, the same tertiary-analysis layer used by clinical genomics teams.

---

## Why this project

This project is divided into two phases. Phase 1 benchmarks GATK and DeepVariant against a real published truth set, measuring precision, recall, and F1 score separately for SNPs and indels. Getting a VCF out of GATK is not enough in clinical genomics, it's important to know how reliable those calls are and, ultimately, what they mean for a patient. Phase 2 takes this further by classifying variants using the actual ACMG/AMP framework, with evidence that can be traced back to its source rather than just assigning a "pathogenic" label.

## Pipeline overview

```
Raw FASTQ (NA12878, chr20)
        │
        ▼
   QC & trimming ── fastp, FastQC, MultiQC
        │
        ▼
   Alignment ── BWA-MEM → samtools sort → GATK MarkDuplicates
        │
        ├──────────────┐
        ▼              ▼
  GATK HaplotypeCaller  DeepVariant
        │              │
        └──────┬───────┘
               ▼
     hap.py benchmark vs. GIAB NA12878 truth set
     (precision / recall / F1, SNP vs. INDEL)
               │
               ▼
     Filter to clinical gene panel (ACMG SF v3.2 default)
               │
               ▼
     VEP annotation
               │
               ▼
     Load into SQLite → query with SQL
               │
               ▼
  ═══════════ PHASE 2 ═══════════
               │
               ▼
     ClinVar + gnomAD + PVS1 evidence extraction
               │
               ▼
     ACMG/AMP 2015 combining-rules classification
               │
               ▼
     Jinja2 HTML clinical report
```

## Tech stack

| Stage | Tools |
|---|---|
| Orchestration | Snakemake |
| Environment | Docker + conda (pinned via `environment.yml`) |
| QC | fastp, FastQC, MultiQC |
| Alignment | BWA-MEM, samtools, GATK MarkDuplicates |
| Variant calling | GATK HaplotypeCaller, Google DeepVariant |
| Benchmarking | hap.py (`xcmp` engine) + GIAB NA12878 truth set |
| Annotation | Ensembl VEP (`--database` mode) |
| Storage/query | SQLite, pysam |
| **Phase 2: Clinical interpretation** | Custom ACMG/AMP 2015 rules engine (Python, pytest-verified) |
| **Phase 2: Evidence sources** | ClinVar, gnomAD (public API, retry/backoff), curated gene-mechanism list |
| **Phase 2: Reporting** | Jinja2 HTML report generator |

## Getting the data

Raw FASTQs aren't committed (too large). Populate `resources/` before running:

```bash
bash scripts/download_reference.sh   # GRCh38 chr20 reference + GIAB NA12878 truth set
bash scripts/build_gene_panel.sh     # ACMG SF v3.2 gene panel BED (edit for a different panel)
```

FASTQ reads: this project uses NA12878 (not the originally-planned HG002 — see "Notes on scope and honesty" below for why), sourced from the GATK public test-data bucket (`gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam`, 5GB, verified real and publicly accessible), subset locally to chr20 with `samtools view`.

## Running the pipeline

```bash
# Build the environment
docker build -t ngs-pipeline .
docker run -v $(pwd):/workspace ngs-pipeline snakemake --cores 4 -p

# Or locally without Docker:
conda env create -f environment.yml
conda env create -f environment-happy.yml   # hap.py needs its own env (old Python 2.7 deps)
conda activate ngs-pipeline
snakemake --cores 4 -p
```

Phase 2, after Phase 1 completes:

```bash
cd clinical-interpretation
pip install pytest requests jinja2
python3 -m pytest tests/ -v                 # verify the classifier against published rules
bash scripts/download_clinvar.sh
python3 germline/classify_variants.py \
    --vcf ../results/annotation/NA12878.panel_filtered.annotated.vcf.gz \
    --clinvar resources/clinvar_chr20.vcf.gz \
    --output results/classifications.csv \
    --query-gnomad
python3 reports/generate_report.py \
    --csv results/classifications.csv \
    --sample-id NA12878 \
    --output reports/NA12878_germline_report.html
```

## Results

### Phase 1: Caller benchmark (GATK vs. DeepVariant, PASS variants, vs. GIAB NA12878 truth)

| Metric | GATK SNP | DeepVariant SNP | GATK INDEL | DeepVariant INDEL |
|---|---|---|---|---|
| Recall | 98.65% | **98.74%** | 94.65% | **96.72%** |
| Precision | 98.76% | **99.31%** | 97.07% | **97.98%** |
| F1 | 98.70% | **99.02%** | 95.84% | **97.35%** |

DeepVariant outperforms GATK across every metric, with the largest gap on indels, consistent with the published literature on DeepVariant's deep-learning approach handling indel-adjacent complexity better than GATK's statistical model. DeepVariant's output was used as the input to Phase 2.

### Phase 2: ACMG/AMP classification (309 variants, 3-gene panel)

| Classification | Count |
|---|---|
| Benign | 201 |
| Uncertain Significance | 108 |

65% of variants resolved to a definitive Benign call using automated ClinVar + gnomAD evidence stacking, up from 0% resolved using ClinVar alone (a single ClinVar-derived Supporting-level code is correctly insufficient on its own under the ACMG/AMP combining rules, population frequency evidence from gnomAD was what actually resolved most of these). No Pathogenic/Likely Pathogenic calls in this healthy reference sample, and zero PVS1 false positives, both expected, correct results for NA12878.

## Customizing the gene panel

The default panel (`resources/panels/acmg_sf_chr20_genes.bed`) targets ACMG SF v3.2 secondary-findings genes on chr20. Swap in a cancer gene panel to reframe this for an oncology-genomics context, update `gene_panel` in `config/config.yaml` and re-run from `filter_gene_panel` onward. Note: Phase 2's PVS1 gene-mechanism list (`clinical-interpretation/germline/pvs1_evidence.py`) would need matching updates for any new gene panel, it's deliberately not automatic.

## Notes on scope and honesty

- Scoped to **chr20 only** and a **3-gene panel** to keep runtime and storage tractable on a laptop — this is a portfolio/demonstration pipeline, not a production clinical pipeline.
- **Sample switched from HG002 to NA12878 mid-project**, after the HG002 300x WGS BAM (126GB) proved unreliable to fetch via remote slicing over an unstable connection — twice. NA12878 is an equally standard GIAB gold-standard sample; the switch is documented here rather than hidden.
- **Population frequency and functional-prediction data (gnomAD, dbNSFP) are queried via public APIs, not downloaded locally**, the full datasets run into the hundreds of GB, wildly disproportionate to a 3-gene panel. This is a deliberate engineering choice, not an oversight; see `clinical-interpretation/README.md` for the full size/access verification behind each data source decision.
- **PVS1 (null-variant) evidence is gated by a small, manually curated gene list** (JAG1, ASXL1), not applied automatically to any null variant, GNAS is deliberately excluded, since its pathogenicity mechanism is imprinting-dependent rather than simple loss-of-function.
- **Automated evidence extraction covers a subset of ACMG/AMP criteria** (PP5, BP6, BA1, BS1, PM2, PVS1). Criteria requiring segregation data, functional studies, or de novo confirmation are not automated.
- The DeepVariant rule shells out to the official Docker image (docker-in-docker).
- Gene panel coordinates were independently verified against Ensembl GRCh38 coordinates after an initial placeholder set was found to contain a GRCh37/GRCh38 mismatch.

## Author

Agnidipa Sett — M.Tech Bioinformatics, Delhi Technological University
[GitHub](https://github.com/Agni-18)
