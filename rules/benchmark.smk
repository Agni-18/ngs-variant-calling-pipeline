SAMPLE_ID = config["sample"]
REF = config["reference"]["fasta"]
TRUTH_VCF = config["truth"]["vcf"]
TRUTH_BED = config["truth"]["bed"]

rule happy_benchmark:
    """
    Benchmark each caller's VCF against the GIAB HG002 truth set
    using hap.py. Produces precision/recall/F1 broken down by
    variant type (SNP vs INDEL) -- this is the honesty check that
    mirrors the LOOCV/DeLong validation approach used elsewhere
    in the portfolio.
    """
    input:
        query=f"results/variants/{SAMPLE_ID}.{{caller}}.vcf.gz",
        tbi=f"results/variants/{SAMPLE_ID}.{{caller}}.vcf.gz.tbi",
        truth=TRUTH_VCF,
        bed=TRUTH_BED,
        ref=REF,
    output:
        summary=f"results/benchmark/{SAMPLE_ID}.{{caller}}.summary.csv",
    params:
        prefix=lambda wc: f"results/benchmark/{SAMPLE_ID}.{wc.caller}",
    log:
        "logs/happy/{caller}.log",
    shell:
        """
        conda run -n happy-env hap.py {input.truth} {input.query} \
            -f {input.bed} \
            -r {input.ref} \
            -o {params.prefix} \
            --engine=xcmp \
            &> {log}
        """


rule benchmark_comparison_plot:
    """
    Combine GATK vs DeepVariant hap.py summaries into a single
    precision/recall/F1 comparison figure for the README.
    """
    input:
        expand("results/benchmark/{sample}.{caller}.summary.csv",
               sample=SAMPLE_ID, caller=["gatk", "deepvariant"]),
    output:
        "results/benchmark/caller_comparison.png",
    script:
        "../scripts/plot_caller_comparison.py"
