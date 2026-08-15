SAMPLE_ID = config["sample"]

rule fastp:
    """Adapter trimming and quality filtering."""
    input:
        r1=config["reads"]["r1"],
        r2=config["reads"]["r2"],
    output:
        r1=f"results/qc/trimmed/{SAMPLE_ID}_R1.trimmed.fastq.gz",
        r2=f"results/qc/trimmed/{SAMPLE_ID}_R2.trimmed.fastq.gz",
        html=f"results/qc/fastp/{SAMPLE_ID}.fastp.html",
        json=f"results/qc/fastp/{SAMPLE_ID}.fastp.json",
    log:
        f"logs/fastp/{SAMPLE_ID}.log",
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.html} -j {output.json} \
              --thread {threads} \
              --detect_adapter_for_pe \
              &> {log}
        """


rule fastqc:
    """Per-sample QC metrics, run post-trim."""
    input:
        r1=f"results/qc/trimmed/{SAMPLE_ID}_R1.trimmed.fastq.gz",
        r2=f"results/qc/trimmed/{SAMPLE_ID}_R2.trimmed.fastq.gz",
    output:
        directory(f"results/qc/fastqc/{SAMPLE_ID}"),
    log:
        f"logs/fastqc/{SAMPLE_ID}.log",
    threads: 2
    shell:
        """
        mkdir -p {output}
        fastqc {input.r1} {input.r2} -o {output} -t {threads} &> {log}
        """


rule multiqc:
    """Aggregate fastp + FastQC into a single report."""
    input:
        f"results/qc/fastqc/{SAMPLE_ID}",
        f"results/qc/fastp/{SAMPLE_ID}.fastp.json",
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc results/qc/ -o results/qc/ -n multiqc_report &> {log}
        """
