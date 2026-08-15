#!/usr/bin/env bash
# Downloads GRCh38 chr20 reference and GIAB HG002 truth set /
# confident regions bed, scoped to chr20 to keep the pipeline
# runnable on a laptop.
set -euo pipefail

mkdir -p resources/reference resources/truth resources/reads resources/panels

echo ">> Downloading GRCh38 chr20 reference..."
curl -L -o resources/reference/chr20.fa.gz \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
gunzip -f resources/reference/chr20.fa.gz
samtools faidx resources/reference/chr20.fa

echo ">> Downloading GIAB HG002 truth VCF + confident regions (chr20)..."
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
curl -L -o resources/truth/HG002_GRCh38_truth.vcf.gz \
    "${GIAB_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
curl -L -o resources/truth/HG002_GRCh38_confident_regions.bed \
    "${GIAB_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

echo ">> Subsetting truth set to chr20 only..."
bcftools view -r chr20 resources/truth/HG002_GRCh38_truth.vcf.gz \
    -Oz -o resources/truth/HG002_GRCh38_chr20_truth.vcf.gz
tabix -p vcf resources/truth/HG002_GRCh38_chr20_truth.vcf.gz
awk '$1 == "chr20"' resources/truth/HG002_GRCh38_confident_regions.bed \
    > resources/truth/HG002_GRCh38_chr20_confident_regions.bed

echo ">> Done. Populate resources/reads/ with HG002 chr20 FASTQs next"
echo "   (see README.md 'Getting the data' section for sources)."
