SAMPLE_ID = config["sample"]
REF = config["reference"]["fasta"]

rule bwa_index:
    """Build BWA index for the reference (one-time, cached)."""
    input:
        REF,
    output:
        f"{REF}.bwt",
        f"{REF}.amb",
        f"{REF}.ann",
        f"{REF}.pac",
        f"{REF}.sa",
    log:
        "logs/bwa_index.log",
    shell:
        "bwa index {input} &> {log}"


rule bwa_mem:
    """Align trimmed reads to the reference, sorted BAM output."""
    input:
        r1=f"results/qc/trimmed/{SAMPLE_ID}_R1.trimmed.fastq.gz",
        r2=f"results/qc/trimmed/{SAMPLE_ID}_R2.trimmed.fastq.gz",
        ref=REF,
        idx=f"{REF}.bwt",
    output:
        bam=f"results/align/{SAMPLE_ID}.sorted.bam",
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA".format(sample=SAMPLE_ID),
    threads: config["threads"]["align"]
    log:
        "logs/bwa_mem.log",
    shell:
        """
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


rule mark_duplicates:
    """Mark PCR/optical duplicates with Picard (or gatk MarkDuplicates)."""
    input:
        bam=f"results/align/{SAMPLE_ID}.sorted.bam",
    output:
        bam=f"results/align/{SAMPLE_ID}.sorted.dedup.bam",
        metrics=f"results/align/{SAMPLE_ID}.dedup_metrics.txt",
    log:
        "logs/mark_duplicates.log",
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            &> {log}
        samtools index {output.bam}
        """


rule alignment_stats:
    """samtools flagstat + coverage summary, feeds into MultiQC."""
    input:
        bam=f"results/align/{SAMPLE_ID}.sorted.dedup.bam",
    output:
        flagstat=f"results/align/{SAMPLE_ID}.flagstat.txt",
        coverage=f"results/align/{SAMPLE_ID}.coverage.txt",
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        samtools coverage {input.bam} > {output.coverage}
        """
