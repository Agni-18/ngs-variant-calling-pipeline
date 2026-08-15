#!/usr/bin/env bash
# Builds a BED file of ACMG SF v3.2 secondary-findings genes that
# fall on chr20, for use as the clinical filtering panel.
# Replace this with a cancer gene panel BED for a MedGenome-oncology
# framing -- e.g. a curated list intersected with a public panel
# such as the TSO500 or a COSMIC Cancer Gene Census subset.
set -euo pipefail

mkdir -p resources/panels

# Example ACMG SF v3.2 genes known to sit on chr20 (illustrative --
# cross-check against the current ACMG SF list before treating this
# as authoritative, the list is versioned and does get updated).
cat > resources/panels/acmg_sf_chr20_genes.bed << 'BED'
chr20	10199206	10236747	JAG1
chr20	32358331	32444075	ASXL1
chr20	57414773	57486247	GNAS
BED

echo ">> Wrote resources/panels/acmg_sf_chr20_genes.bed"
echo "   NOTE: coordinates are illustrative placeholders -- pull exact"
echo "   GRCh38 coordinates from Ensembl/UCSC before running the real"
echo "   pipeline. See README 'Customizing the gene panel'."
