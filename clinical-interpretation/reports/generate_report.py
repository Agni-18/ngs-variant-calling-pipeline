"""
Generate an HTML germline variant classification report from the
CSV output of classify_variants.py.
"""
import argparse
import csv
from collections import Counter
from datetime import date
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

CLASSIFICATION_SLUGS = {
    "Pathogenic": "pathogenic",
    "Likely Pathogenic": "likely-pathogenic",
    "Uncertain Significance": "uncertain-significance",
    "Likely Benign": "likely-benign",
    "Benign": "benign",
}


def load_variants(csv_path: str) -> list:
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        variants = list(reader)

    for v in variants:
        v["classification_slug"] = CLASSIFICATION_SLUGS.get(
            v["classification"], "uncertain-significance"
        )
        if v.get("gnomad_af"):
            try:
                af = float(v["gnomad_af"])
                v["gnomad_af_display"] = f"{af:.4f}" if af >= 0.0001 else f"{af:.2e}"
            except ValueError:
                v["gnomad_af_display"] = v["gnomad_af"]
        else:
            v["gnomad_af_display"] = "\u2014"

    return variants


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference-genome", default="GRCh38")
    parser.add_argument("--region", default="chr20")
    parser.add_argument("--gene-panel", default="JAG1, ASXL1, GNAS (ACMG SF v3.2 subset)")
    parser.add_argument("--variant-caller", default="DeepVariant 1.6.0")
    args = parser.parse_args()

    variants = load_variants(args.csv)
    counts = Counter(v["classification"] for v in variants)

    priority = {
        "Pathogenic": 0, "Likely Pathogenic": 1,
        "Uncertain Significance": 2,
        "Likely Benign": 3, "Benign": 4,
    }
    variants.sort(key=lambda v: priority.get(v["classification"], 5))

    template_dir = Path(__file__).resolve().parents[1] / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    template = env.get_template("germline_report.html")

    html = template.render(
        sample_id=args.sample_id,
        reference_genome=args.reference_genome,
        region=args.region,
        gene_panel=args.gene_panel,
        variant_caller=args.variant_caller,
        generated_date=date.today().isoformat(),
        counts=counts,
        variants=variants,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)

    print(f"Report written to {output_path}")
    print(f"  {len(variants)} variants, breakdown: {dict(counts)}")


if __name__ == "__main__":
    main()
