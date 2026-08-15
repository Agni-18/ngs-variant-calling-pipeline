"""
NGS Variant Calling & Clinical Annotation Pipeline
===================================================
Sample: GIAB HG002 (chr20 subset)
Author: Agnidipa Sett

Pipeline stages:
  1. QC            -> fastp, FastQC, MultiQC
  2. Alignment      -> BWA-MEM, samtools sort/dedup
  3. Variant calling -> GATK HaplotypeCaller AND DeepVariant (compare both)
  4. Benchmarking   -> hap.py against GIAB truth VCF (precision/recall/F1)
  5. Annotation     -> VEP, filtered to a clinical gene panel
  6. Storage        -> load annotated variants into SQLite for querying

Run:
    docker build -t ngs-pipeline .
    docker run -v $(pwd):/workspace ngs-pipeline snakemake --cores 4 --use-conda
"""

configfile: "config/config.yaml"

SAMPLE = config["sample"]
CHROM = config["chromosome"]
CALLERS = ["gatk", "deepvariant"]

include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/variant_calling.smk"
include: "rules/benchmark.smk"
include: "rules/annotate.smk"
include: "rules/database.smk"


rule all:
    input:
        # QC reports
        f"results/qc/multiqc_report.html",
        # Alignment
        f"results/align/{SAMPLE}.sorted.dedup.bam",
        # Variant calls from both callers
        expand("results/variants/{sample}.{caller}.vcf.gz", sample=SAMPLE, caller=CALLERS),
        # Benchmark summary comparing both callers against GIAB truth
        expand("results/benchmark/{sample}.{caller}.summary.csv", sample=SAMPLE, caller=CALLERS),
        # Annotated, panel-filtered variants
        f"results/annotation/{SAMPLE}.panel_filtered.annotated.vcf.gz",
        # SQLite database of final variants
        f"results/db/{SAMPLE}_variants.sqlite"


rule clean:
    shell:
        "rm -rf results/*"
