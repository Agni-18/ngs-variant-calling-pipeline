# Clinical Variant Calling Pipeline — GIAB HG002 (chr20)

A reproducible, benchmarked WGS/WES variant calling and annotation pipeline built with **Snakemake** and **Docker**, validated against the **Genome in a Bottle (GIAB) HG002** truth set.

Two variant callers — **GATK HaplotypeCaller** and **Google DeepVariant** — are run in parallel and benchmarked head-to-head with **hap.py**, rather than trusting a single caller's output. Final variants are annotated with **VEP**, filtered to a clinical gene panel, and loaded into **SQLite** for querying.

> Built as a portfolio project for bioinformatics/computational genomics roles (clinical NGS pipelines, variant interpretation). Scope is deliberately narrowed to chromosome 20 so the full pipeline — including truth-set benchmarking — runs on a laptop.

---

## Why this project

Most student variant-calling pipelines stop at "I ran GATK and got a VCF." This one asks the next question: **how good are those calls, actually?** Benchmarking against a published truth set with precision/recall/F1 broken out by SNP vs. indel is the same kind of honest validation used in [my T2D transcriptomics signature work](#) (leakage-corrected LOOCV, DeLong's test, permutation testing) — the goal is calls you can defend, not just calls you can produce.

## Pipeline overview

```
Raw FASTQ (HG002, chr20)
        │
        ▼
   QC & trimming ── fastp, FastQC, MultiQC
        │
        ▼
   Alignment ── BWA-MEM → samtools sort → Picard/GATK MarkDuplicates
        │
        ├──────────────┐
        ▼              ▼
  GATK HaplotypeCaller  DeepVariant
        │              │
        └──────┬───────┘
               ▼
     hap.py benchmark vs. GIAB HG002 truth set
     (precision / recall / F1, SNP vs. INDEL)
               │
               ▼
     Best caller's VCF → VEP annotation
               │
               ▼
     Filter to clinical gene panel (ACMG SF v3.2 default)
               │
               ▼
     Load into SQLite → query with SQL
```

## Tech stack

| Stage | Tools |
|---|---|
| Orchestration | Snakemake |
| Environment | Docker + conda (pinned via `environment.yml`) |
| QC | fastp, FastQC, MultiQC |
| Alignment | BWA-MEM, samtools, GATK MarkDuplicates |
| Variant calling | GATK HaplotypeCaller, Google DeepVariant |
| Benchmarking | hap.py + GIAB HG002 truth set (vcfeval engine) |
| Annotation | Ensembl VEP |
| Storage/query | SQLite, pysam |
| CI | GitHub Actions (Snakemake lint + dry-run) |

## Repository structure

```
.
├── Snakefile                 # Pipeline entrypoint (rule all)
├── config/config.yaml         # Sample, reference, panel, thread config
├── rules/
│   ├── qc.smk
│   ├── align.smk
│   ├── variant_calling.smk
│   ├── benchmark.smk
│   ├── annotate.smk
│   └── database.smk
├── scripts/
│   ├── download_reference.sh  # Pulls chr20 ref + GIAB truth set
│   ├── build_gene_panel.sh    # Builds clinical filter panel BED
│   ├── load_to_sqlite.py      # VCF -> SQLite ETL
│   └── plot_caller_comparison.py
├── Dockerfile
├── environment.yml
└── .github/workflows/lint.yml
```

## Getting the data

Raw FASTQs aren't committed to the repo (too large). Populate `resources/` before running:

```bash
# Reference genome + GIAB truth set, subset to chr20
bash scripts/download_reference.sh

# Clinical gene panel (edit scripts/build_gene_panel.sh for a
# disease-specific panel, e.g. an oncology panel)
bash scripts/build_gene_panel.sh

# HG002 chr20 FASTQs: subset from the full GIAB HG002 Illumina
# WGS reads (NIST/GIAB FTP or the precisionFDA HG002 dataset),
# e.g. using `samtools view -b <bam> chr20` on the public HG002
# BAM, then `samtools fastq` to regenerate paired FASTQs.
# Place the result at resources/reads/HG002_chr20_R{1,2}.fastq.gz
```

## Running the pipeline

```bash
# Build the environment
docker build -t ngs-pipeline .

# Dry run (check the DAG without executing)
docker run -v $(pwd):/workspace ngs-pipeline snakemake -n

# Full run
docker run -v $(pwd):/workspace ngs-pipeline snakemake --cores 4 -p

# Or locally without Docker, using conda directly:
conda env create -f environment.yml
conda activate ngs-pipeline
snakemake --cores 4 -p
```

## Results

*(Fill in after running on real HG002 chr20 data — this is the section that shows benchmarking rigor, keep it here rather than burying it in results/.)*

- **Caller comparison (precision / recall / F1, PASS variants):**
  `results/benchmark/caller_comparison.png`
- **Chosen caller for downstream annotation:** *TBD after inspecting hap.py output*
- **Annotated, panel-filtered variant count:** *TBD*
- **Example SQL query** (pathogenic-flagged variants in panel genes):

```sql
SELECT chrom, pos, ref, alt, gene_symbol, consequence, clin_sig
FROM variants
WHERE clin_sig LIKE '%pathogenic%'
ORDER BY chrom, pos;
```

## Customizing the gene panel

The default panel (`resources/panels/acmg_sf_chr20_genes.bed`) targets ACMG SF v3.2 secondary-findings genes on chr20. Swap in a cancer gene panel (e.g. a Cancer Gene Census subset, or a clinical oncology panel like TSO500) to reframe this pipeline for an oncology-genomics context — update `gene_panel` in `config/config.yaml` and re-run from the `filter_gene_panel` rule onward.

## Notes on scope and honesty

- Scoped to **chr20 only** to keep runtime and storage tractable on a laptop — this is a portfolio/demonstration pipeline, not a production clinical pipeline (no multi-sample joint genotyping, no CNV/SV calling, no full-genome QC thresholds).
- The DeepVariant rule shells out to the official Docker image (docker-in-docker); swap for the bioconda `deepvariant` package if you'd rather avoid nested Docker.
- Gene panel coordinates in `build_gene_panel.sh` are illustrative placeholders — verify against current Ensembl/UCSC GRCh38 coordinates and the current ACMG SF list before treating results as clinically meaningful.

## Author

Agnidipa Sett — M.Tech Bioinformatics, Delhi Technological University
[GitHub](https://github.com/Agni-18)
