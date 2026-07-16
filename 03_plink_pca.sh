#!/usr/bin/env bash
# =============================================================================
# 03_plink_pca.sh
# Runs PLINK PCA for all four MAF-filtered single-SNP datasets.
#
# LD pruning is NOT performed — one SNP per RAD locus was already selected
# at the populations stage (--write-random-snp), which provides sufficient
# LD control for 2bRAD data. Between-locus LD is expected to be minimal.
#
# MAF filter (>0.05) is applied in PLINK for A2-A4 datasets (within-taxon
# runs that were exported without MAF filtering from populations, to preserve
# flexibility). A1 was already MAF-filtered at the populations stage but we
# apply the PLINK MAF filter consistently across all four datasets.
#
# Requirements (conda env: pca):
#   plink, bcftools, htslib (bgzip/tabix)
#
# Run from project root.
# =============================================================================
set -euo pipefail

POPBASE="02_datasets/populations_runs"
OUTBASE="03_analyses/pca"
LOGDIR="03_analyses/pca/logs"
THREADS=4

mkdir -p "${LOGDIR}"

# ── datasets to process ───────────────────────────────────────────────────
# Format: "tag|vcf_subpath"
# VCF path is relative to POPBASE/<tag>/
DATASETS=(
    "A1_all_taxa_MAF_randomSNP|populations.snps.fixed.vcf.gz"
    "A2_Ss_F1F2_noMAF_randomSNP|populations.snps.fixed.vcf.gz"
    "A3_Si_wildstock_noMAF_randomSNP|populations.snps.fixed.vcf.gz"
    "A4_Sk_F1F2_noMAF_randomSNP|populations.snps.fixed.vcf.gz"
)

# ── helper ────────────────────────────────────────────────────────────────
run_pca() {
    local tag="$1"
    local vcf="$2"

    local outdir="${OUTBASE}/${tag}"
    local logfile="${LOGDIR}/plink_${tag}.log"
    mkdir -p "${outdir}"

    local tmpdir
    tmpdir="$(mktemp -d)"
    trap 'rm -rf "${tmpdir}"' RETURN

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  TAG : ${tag}"
    echo "  VCF : ${vcf}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # ── Step 1: rename numeric CHROM values to cN format ─────────────────
    # PLINK requires chromosome names that look like chromosomes.
    # Stacks outputs numeric contig IDs (e.g. 1, 29, 347); rename to c1, c29, c347.
    local rename_map="${tmpdir}/rename_chrs.txt"
    local vcf_renamed="${tmpdir}/renamed.vcf.gz"

    bcftools query -f '%CHROM\n' "${vcf}" | sort -u \
        | awk '{print $1 "\tc" $1}' > "${rename_map}"

    bcftools annotate \
        --rename-chrs "${rename_map}" \
        -Oz -o "${vcf_renamed}" \
        "${vcf}"
    tabix -f -p vcf "${vcf_renamed}"

    # ── Step 2: VCF → PLINK bed ───────────────────────────────────────────
    # --double-id: use sample ID as both FID and IID (no family structure)
    # --allow-extra-chr: permit non-standard chromosome names (cN format)
    # --set-missing-var-ids: assign variant IDs from CHROM:POS:REF:ALT
    # --maf 0.05: MAF filter applied consistently across all datasets
    #             (A1 was pre-filtered in populations; A2-A4 were not)
    plink \
        --vcf "${vcf_renamed}" \
        --double-id \
        --allow-extra-chr \
        --set-missing-var-ids '@:#:$1:$2' \
        --maf 0.05 \
        --make-bed \
        --out "${outdir}/plink_input" \
        2>&1 | tee "${logfile}"

    # ── Step 3: PCA ───────────────────────────────────────────────────────
    # Run 20 PCs so PC axis arguments in the R script can go up to PC20.
    plink \
        --bfile "${outdir}/plink_input" \
        --allow-extra-chr \
        --pca 20 \
        --out "${outdir}/pca" \
        2>&1 | tee -a "${logfile}"

    echo "  PCA complete → ${outdir}/"
    echo "    eigenvec : ${outdir}/pca.eigenvec"
    echo "    eigenval : ${outdir}/pca.eigenval"
    echo ""
}

# ── main loop ─────────────────────────────────────────────────────────────
for entry in "${DATASETS[@]}"; do
    tag="${entry%%|*}"
    vcf_rel="${entry##*|}"
    vcf="${POPBASE}/${tag}/${vcf_rel}"

    if [[ ! -s "${vcf}" ]]; then
        echo "SKIP (VCF not found): ${vcf}"
        continue
    fi

    run_pca "${tag}" "${vcf}"
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "All PCA runs complete. Outputs in: ${OUTBASE}/"
echo ""
echo "Next step: Rscript scripts/03_plot_pca.R"
