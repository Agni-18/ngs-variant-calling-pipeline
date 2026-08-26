#!/usr/bin/env bash
set -euo pipefail
mkdir -p resources
echo ">> Downloading ClinVar GRCh38 VCF (~184MB)..."
curl -L -o resources/clinvar.vcf.gz \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
echo ">> Indexing..."
tabix -p vcf resources/clinvar.vcf.gz
echo ">> Subsetting to chr20 (ClinVar uses '20', not 'chr20')..."
bcftools view -r 20 resources/clinvar.vcf.gz -Oz -o resources/clinvar_chr20.vcf.gz
tabix -p vcf resources/clinvar_chr20.vcf.gz
echo ">> Done. Variant count:"
zcat resources/clinvar_chr20.vcf.gz | grep -vc "^#"
