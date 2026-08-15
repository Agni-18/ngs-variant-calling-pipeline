"""
Parse an annotated, panel-filtered VCF and load it into a SQLite
database as a tidy variants table.

Called as a Snakemake `script:` directive, so `snakemake` object
(input/output/log) is injected automatically -- see rules/database.smk.
"""
import sqlite3
import sys

import pysam

vcf_path = snakemake.input.vcf
db_path = snakemake.output.db
log_path = snakemake.log[0]

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS variants (
    chrom TEXT,
    pos INTEGER,
    ref TEXT,
    alt TEXT,
    qual REAL,
    filter TEXT,
    gene_symbol TEXT,
    consequence TEXT,
    impact TEXT,
    clin_sig TEXT,
    gnomad_af REAL,
    genotype TEXT
);
"""


def extract_vep_field(info_csq, field, csq_format):
    """Pull a single field (e.g. 'SYMBOL', 'IMPACT') out of VEP's CSQ string."""
    if info_csq is None:
        return None
    fields = csq_format.split("|")
    if field not in fields:
        return None
    idx = fields.index(field)
    # CSQ can carry multiple transcript annotations, comma-separated;
    # take the first (VEP typically orders by severity with --pick,
    # otherwise this should be revisited depending on annotation flags used)
    first_csq = info_csq.split(",")[0]
    parts = first_csq.split("|")
    return parts[idx] if idx < len(parts) else None


def main():
    with open(log_path, "w") as log:
        vcf = pysam.VariantFile(vcf_path)

        # Extract the CSQ format string VEP writes into the VCF header
        csq_format = ""
        for record in vcf.header.records:
            if record.key == "INFO" and record.get("ID") == "CSQ":
                desc = record.get("Description", "")
                if "Format:" in desc:
                    csq_format = desc.split("Format:")[1].strip().strip('"')

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute(CREATE_TABLE_SQL)

        n_variants = 0
        for rec in vcf.fetch():
            csq = rec.info.get("CSQ")
            csq_str = csq[0] if isinstance(csq, tuple) else csq

            gene = extract_vep_field(csq_str, "SYMBOL", csq_format)
            consequence = extract_vep_field(csq_str, "Consequence", csq_format)
            impact = extract_vep_field(csq_str, "IMPACT", csq_format)
            clin_sig = extract_vep_field(csq_str, "CLIN_SIG", csq_format)
            gnomad_af = extract_vep_field(csq_str, "gnomAD_AF", csq_format)

            for sample_name in rec.samples:
                gt = rec.samples[sample_name]["GT"]
                genotype = "/".join(str(a) for a in gt) if gt else None

                cur.execute(
                    """
                    INSERT INTO variants
                    (chrom, pos, ref, alt, qual, filter, gene_symbol,
                     consequence, impact, clin_sig, gnomad_af, genotype)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        rec.chrom,
                        rec.pos,
                        rec.ref,
                        ",".join(rec.alts) if rec.alts else None,
                        rec.qual,
                        ";".join(rec.filter.keys()) if rec.filter else None,
                        gene,
                        consequence,
                        impact,
                        clin_sig,
                        float(gnomad_af) if gnomad_af not in (None, "") else None,
                        genotype,
                    ),
                )
            n_variants += 1

        conn.commit()
        conn.close()
        log.write(f"Loaded {n_variants} variants into {db_path}\n")


if __name__ == "__main__":
    main()
