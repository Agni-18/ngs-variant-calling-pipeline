"""
Run ACMG/AMP classification on the Phase 1 pipeline's annotated,
panel-filtered VCF, using ClinVar, gnomAD, and gene-gated PVS1
evidence.

Usage:
    python3 classify_variants.py \
        --vcf ../results/annotation/NA12878.panel_filtered.annotated.vcf.gz \
        --clinvar resources/clinvar_chr20.vcf.gz \
        --output results/classifications.csv \
        --query-gnomad
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from germline.classifier import classify
from germline.clinvar_evidence import get_clinvar_evidence, load_clinvar_index
from germline.gnomad_evidence import get_frequency_evidence, query_gnomad_af
from germline.pvs1_evidence import get_pvs1_evidence


def parse_csq_format(vcf_path: str) -> list[str]:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                if "Format:" in line:
                    format_str = line.split("Format:")[1].strip().strip('">')
                    return format_str.split("|")
            if not line.startswith("#"):
                break
    return []


def extract_gene_symbol(info_dict: dict, csq_fields: list[str]) -> str:
    if "CSQ" not in info_dict or "SYMBOL" not in csq_fields:
        return ""
    symbol_idx = csq_fields.index("SYMBOL")
    first_csq = info_dict["CSQ"].split(",")[0]
    parts = first_csq.split("|")
    return parts[symbol_idx] if symbol_idx < len(parts) else ""


def extract_consequence(info_dict: dict, csq_fields: list[str]) -> str:
    if "CSQ" not in info_dict or "Consequence" not in csq_fields:
        return ""
    consequence_idx = csq_fields.index("Consequence")
    first_csq = info_dict["CSQ"].split(",")[0]
    parts = first_csq.split("|")
    return parts[consequence_idx] if consequence_idx < len(parts) else ""


def parse_vcf_variants(vcf_path: str):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            info_dict = {}
            for entry in fields[7].split(";"):
                if "=" in entry:
                    key, _, value = entry.partition("=")
                    info_dict[key] = value
            yield chrom, pos, ref, alt, info_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Phase 1 annotated VCF")
    parser.add_argument("--clinvar", required=True, help="chr20-subset ClinVar VCF")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--query-gnomad", action="store_true",
        help="Also query gnomAD's public API for population allele frequency evidence.")
    parser.add_argument("--limit", type=int, default=None,
        help="Process only the first N variant records (for testing).")
    args = parser.parse_args()

    print(f"Loading ClinVar index from {args.clinvar}...", file=sys.stderr)
    clinvar_index = load_clinvar_index(args.clinvar)
    print(f"  {len(clinvar_index)} ClinVar records loaded", file=sys.stderr)

    csq_fields = parse_csq_format(args.vcf)
    print(f"CSQ field order detected: {csq_fields[:5]}...", file=sys.stderr)

    rows = []
    n_variants = 0
    n_in_clinvar = 0

    for chrom, pos, ref, alt, info_dict in parse_vcf_variants(args.vcf):
        if args.limit is not None and n_variants >= args.limit:
            break
        n_variants += 1

        for single_alt in alt.split(","):
            gene_symbol = extract_gene_symbol(info_dict, csq_fields)
            consequence = extract_consequence(info_dict, csq_fields)

            evidence_codes, clinvar_record = get_clinvar_evidence(
                chrom, pos, ref, single_alt, clinvar_index
            )
            if clinvar_record is not None:
                n_in_clinvar += 1

            evidence_codes = evidence_codes + get_pvs1_evidence(gene_symbol, consequence)

            gnomad_af = None
            if args.query_gnomad:
                gnomad_af = query_gnomad_af(chrom, pos, ref, single_alt)
                evidence_codes = evidence_codes + get_frequency_evidence(gnomad_af)

            result = classify(evidence_codes)

            rows.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": single_alt,
                "gene_symbol": gene_symbol,
                "consequence": consequence,
                "evidence_codes": ";".join(evidence_codes) if evidence_codes else "(none)",
                "classification": result.classification.value,
                "rule_matched": result.rule_matched,
                "gnomad_af": gnomad_af if gnomad_af is not None else "",
                "clinvar_clnsig": clinvar_record.clnsig if clinvar_record else "",
                "clinvar_review_status": clinvar_record.clnrevstat if clinvar_record else "",
                "clinvar_disease": clinvar_record.clndn if clinvar_record else "",
            })

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nProcessed {n_variants} variant records.", file=sys.stderr)
    print(f"{n_in_clinvar} alleles matched an existing ClinVar entry.", file=sys.stderr)
    print(f"Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
