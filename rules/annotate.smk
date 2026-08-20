SAMPLE_ID = config["sample"]
GENE_PANEL = config["gene_panel"]
ASSEMBLY = config["annotation"]["assembly"]

BEST_CALLER = "deepvariant"

rule filter_gene_panel:
    input:
        vcf=f"results/variants/{SAMPLE_ID}.{BEST_CALLER}.vcf.gz",
        tbi=f"results/variants/{SAMPLE_ID}.{BEST_CALLER}.vcf.gz.tbi",
        panel=GENE_PANEL,
    output:
        vcf=f"results/variants/{SAMPLE_ID}.panel_filtered.vcf.gz",
    log:
        "logs/filter_gene_panel.log",
    shell:
        """
        bcftools view -R {input.panel} {input.vcf} -Oz -o {output.vcf} &> {log}
        tabix -p vcf {output.vcf}
        """




rule vep_annotate:
    input:
        vcf=f"results/variants/{SAMPLE_ID}.panel_filtered.vcf.gz",
    output:
        vcf=f"results/annotation/{SAMPLE_ID}.panel_filtered.annotated.vcf.gz",
        html=f"results/annotation/{SAMPLE_ID}.panel_filtered.annotated.vcf.gz_summary.html",
    threads: config["threads"]["annotate"]
    log:
        "logs/vep.log",
    shell:
        """
        vep --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf --compress_output bgzip \
            --database \
            --assembly {ASSEMBLY} \
            --symbol --biotype --canonical \
            --hgvs --check_existing \
            --numbers --variant_class \
            --force_overwrite \
            &> {log}
        tabix -p vcf {output.vcf}
        """
