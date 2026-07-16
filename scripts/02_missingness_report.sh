#!/usr/bin/env bash
# =============================================================================
# 02_missingness_report.sh
# Generate per-individual and per-site missingness reports for key datasets.
# Review these before running any analyses to catch problematic samples.
#
# Requirements (conda env: vcftools):
#   vcftools
#
# Run from project root.
# =============================================================================
set -euo pipefail

POPBASE="02_datasets/populations_runs"
REPDIR="02_datasets/missingness_reports"
mkdir -p "${REPDIR}"

report_missingness() {
    local tag="$1"
    local vcf="${POPBASE}/${tag}/populations.snps.vcf"
    local outprefix="${REPDIR}/${tag}"

    if [[ ! -s "${vcf}" ]]; then
        echo "  SKIP (no VCF): ${tag}"
        return
    fi

    echo "  ${tag}"
    vcftools --vcf "${vcf}" --missing-indv --out "${outprefix}" 2>/dev/null
    vcftools --vcf "${vcf}" --missing-site  --out "${outprefix}" 2>/dev/null
    echo "    → ${outprefix}.imiss  (per-individual)"
    echo "    → ${outprefix}.lmiss  (per-site)"
}

# Run reports for all key datasets
# Focus on the between-taxa and within-taxon single-SNP sets
# (the allSNP sets will have higher missingness by design)
TAGS=(
    "A1_all_taxa_MAF_randomSNP"
    "A1_all_taxa_noMAF_randomSNP"
    "A2_Ss_F1F2_noMAF_randomSNP"
    "A3_Si_wildstock_noMAF_randomSNP"
    "A4_Sk_F1F2_noMAF_randomSNP"
)

echo "Generating missingness reports..."
for tag in "${TAGS[@]}"; do
    report_missingness "${tag}"
done

echo ""
echo "Reports written to: ${REPDIR}/"
echo ""
echo "Review .imiss files for samples with F_MISS > 0.30."
echo "Flag any outlier samples before proceeding to analyses."
echo "Remember to check Si_HAPAMA14_04 especially (generation=unknown→stock)."
