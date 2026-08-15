SAMPLE_ID = config["sample"]
REF = config["reference"]["fasta"]
CHROM = config["chromosome"]

rule gatk_haplotypecaller:
    """Germline short variant calling with GATK HaplotypeCaller."""
    input:
        bam=f"results/align/{SAMPLE_ID}.sorted.dedup.bam",
        ref=REF,
    output:
        vcf=f"results/variants/{SAMPLE_ID}.gatk.vcf.gz",
    threads: config["threads"]["call"]
    log:
        "logs/gatk_haplotypecaller.log",
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {CHROM} \
            -O {output.vcf} \
            &> {log}
        """


rule deepvariant:
    """
    Germline short variant calling with Google DeepVariant.
    Run via the official Docker image so no local install is needed
    beyond Docker itself (see Dockerfile for the base pipeline image
    which already has docker-in-docker / singularity configured).
    """
    input:
        bam=f"results/align/{SAMPLE_ID}.sorted.dedup.bam",
        ref=REF,
    output:
        vcf=f"results/variants/{SAMPLE_ID}.deepvariant.vcf.gz",
    threads: config["threads"]["call"]
    log:
        "logs/deepvariant.log",
    shell:
        """
        docker run --rm \
            -v "$(pwd)":/workspace \
            google/deepvariant:1.6.0 \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=/workspace/{input.ref} \
            --reads=/workspace/{input.bam} \
            --regions={CHROM} \
            --output_vcf=/workspace/{output.vcf} \
            --num_shards={threads} \
            &> {log}
        """


rule sort_index_vcf:
    """Ensure both caller outputs are bgzipped + tabix-indexed for hap.py."""
    input:
        vcf="results/variants/{sample}.{caller}.vcf.gz",
    output:
        tbi="results/variants/{sample}.{caller}.vcf.gz.tbi",
    shell:
        "tabix -p vcf {input.vcf}"
