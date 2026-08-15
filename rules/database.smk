SAMPLE_ID = config["sample"]

rule load_to_sqlite:
    """
    Parse the final annotated, panel-filtered VCF into a tidy
    SQLite table so variants can be explored with plain SQL --
    useful both as a portfolio SQL demo and as a lightweight
    stand-in for a clinical variant database.
    """
    input:
        vcf=f"results/annotation/{SAMPLE_ID}.panel_filtered.annotated.vcf.gz",
    output:
        db=config["database"]["path"],
    log:
        "logs/load_to_sqlite.log",
    script:
        "../scripts/load_to_sqlite.py"
