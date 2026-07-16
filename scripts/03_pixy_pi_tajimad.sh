#!/usr/bin/env bash
# =============================================================================
# 03_pixy_pi_tajimad.sh
# Runs pixy to compute Pi, Dxy, Fst, and Tajima's D across all datasets.
#
# ANALYSES:
#   Pi / Dxy / Fst:
#     A1: Ss vs Sk vs Si        (popmap_all_taxa_species.tsv)
#     A2: Ss PAKA(F1) vs PAHY(F2) (popmap_Ss_F1F2.tsv)
#     A3: Si wild vs stock      (popmap_Si_wildstock.tsv)
#     A4: Sk F1 vs F2           (popmap_Sk_F1F2.tsv)
#
#   Tajima's D (all-taxa VCF only):
#     Run 1: Ss / Sk / Si separate  (popmap_all_taxa_species.tsv)
#     Run 2: Ss+Sk pooled as SsSk   (popmap_SsSk_pooled_Si.tsv)
#            Tests Wahlund-SFS effect -- pooled D should be more positive
#            than either taxon separately if they are genuinely diverged.
#
# WINDOW SIZE: 36 bp (= one 2bRAD locus per window; Stacks contigs are
#   single loci so each contig contains one RAD tag)
#
# INPUT VCFs: IUPAC-coded with invariant sites from phylip-to-vcf converter.
#   This script bgzips and tabix-indexes them before running pixy.
#
# Requirements (conda env: pixy): pixy, bgzip, tabix
# Run from project root.
# =============================================================================
set -euo pipefail

PHYLIP_FIXED="02_datasets/phylip_fixed"
POPMAP_DIR="metadata/popmaps"
OUTBASE="03_analyses/pixy"
LOGDIR="03_analyses/pixy/logs"
THREADS=4
WINDOW=36

mkdir -p "${OUTBASE}" "${LOGDIR}"

# =============================================================================
# STEP 1 — bgzip + tabix index all invariant VCFs
# =============================================================================

# echo "================================================================"
# echo "STEP 1: Compressing and indexing invariant VCFs"
# echo "================================================================"

# index_vcf() {
#     local vcf="$1"
#     local vcf_gz="${vcf}.gz"

#     if [[ ! -s "${vcf}" ]]; then
#         echo "  ERROR: VCF not found: ${vcf}"
#         echo "  Ensure phylip-to-vcf conversion has been run first."
#         exit 1
#     fi

#     if [[ -s "${vcf_gz}" ]]; then
#         echo "  Already compressed, re-indexing: ${vcf_gz}"
#         tabix -f -p vcf "${vcf_gz}"
#     else
#         echo "  Compressing: ${vcf}"
#         bgzip -f "${vcf}"
#         tabix -f -p vcf "${vcf_gz}"
#     fi
#     echo "    -> ${vcf_gz}"
# }

VCF_A1="${PHYLIP_FIXED}/A1_all_taxa_indiv_phylip/A1_all_taxa_indiv_invariants.vcf"
VCF_A2="${PHYLIP_FIXED}/A2_Ss_indiv_phylip/A2_Ss_indiv_invariants.vcf"
VCF_A3="${PHYLIP_FIXED}/A3_Si_indiv_phylip/A3_Si_indiv_invariants.vcf"
VCF_A4="${PHYLIP_FIXED}/A4_Sk_indiv_phylip/A4_Sk_indiv_invariants.vcf"

# index_vcf "${VCF_A1}"
# index_vcf "${VCF_A2}"
# index_vcf "${VCF_A3}"
# index_vcf "${VCF_A4}"

# After bgzip the .vcf.gz paths are used for all pixy runs
VCF_A1_GZ="${VCF_A1}.gz"
VCF_A2_GZ="${VCF_A2}.gz"
VCF_A3_GZ="${VCF_A3}.gz"
VCF_A4_GZ="${VCF_A4}.gz"

echo ""

# =============================================================================
# STEP 2 — pixy helper
# =============================================================================

run_pixy() {
    local run_tag="$1"
    local vcf_gz="$2"
    local popmap="$3"
    local stats="$4"
    local outdir="${OUTBASE}/${run_tag}"
    local logfile="${LOGDIR}/pixy_${run_tag}.log"

    mkdir -p "${outdir}"

    echo "----------------------------------------------------------------"
    echo "  Run    : ${run_tag}"
    echo "  Stats  : ${stats}"
    echo "  Popmap : ${popmap}"
    echo "  Window : ${WINDOW} bp"

    if [[ ! -s "${vcf_gz}" ]]; then
        echo "  ERROR: VCF not found: ${vcf_gz}"; return 1
    fi
    if [[ ! -s "${popmap}" ]]; then
        echo "  ERROR: Popmap not found: ${popmap}"; return 1
    fi

    pixy \
        --stats "${stats}" \
        --vcf "${vcf_gz}" \
        --populations "${popmap}" \
        --window_size "${WINDOW}" \
        --n_cores "${THREADS}" \
        --output_folder "${outdir}" \
        --output_prefix "${run_tag}" \
        2>&1 | tee "${logfile}"

    echo "  Done -> ${outdir}/"
    echo ""
}

# =============================================================================
# STEP 3 — Pi / Dxy / Fst runs
# =============================================================================

echo "================================================================"
echo "STEP 2: Pi, Dxy, Fst"
echo "================================================================"
echo ""

run_pixy \
    "A1_all_taxa_pi" \
    "${VCF_A1_GZ}" \
    "${POPMAP_DIR}/popmap_all_taxa_species.tsv" \
    'pi' 'dxy' 'fst'

run_pixy \
    "A2_Ss_F1F2_pi" \
    "${VCF_A2_GZ}" \
    "${POPMAP_DIR}/popmap_Ss_F1F2.tsv" \
    'pi' 'dxy' 'fst'

run_pixy \
    "A3_Si_wildstock_pi" \
    "${VCF_A3_GZ}" \
    "${POPMAP_DIR}/popmap_Si_wildstock.tsv" \
    'pi' 'dxy' 'fst'

run_pixy \
    "A4_Sk_F1F2_pi" \
    "${VCF_A4_GZ}" \
    "${POPMAP_DIR}/popmap_Sk_F1F2.tsv" \
    'pi' 'dxy' 'fst'

# =============================================================================
# STEP 4 — Tajima's D runs
# =============================================================================

echo "================================================================"
echo "STEP 3: Tajima's D (all-taxa VCF)"
echo "================================================================"
echo ""

# Run 1: Ss / Sk / Si as separate populations
run_pixy \
    "A1_all_taxa_tajimad_separate" \
    "${VCF_A1_GZ}" \
    "${POPMAP_DIR}/popmap_all_taxa_species.tsv" \
    'tajima_d'

# Run 2: Ss+Sk pooled as SsSk, Si separate
# A positive shift in pooled D vs. separate Ss/Sk D reflects the
# Wahlund-SFS effect expected when two diverged gene pools are merged
run_pixy \
    "A1_SsSk_pooled_tajimad" \
    "${VCF_A1_GZ}" \
    "${POPMAP_DIR}/popmap_SsSk_pooled_Si.tsv" \
    'tajima_d'

# =============================================================================
# STEP 5 — Summary
# =============================================================================

echo "================================================================"
echo "All pixy runs complete."
echo ""
echo "Output directories:"
echo "  ${OUTBASE}/A1_all_taxa_pi/"
echo "  ${OUTBASE}/A2_Ss_F1F2_pi/"
echo "  ${OUTBASE}/A3_Si_wildstock_pi/"
echo "  ${OUTBASE}/A4_Sk_F1F2_pi/"
echo "  ${OUTBASE}/A1_all_taxa_tajimad_separate/"
echo "  ${OUTBASE}/A1_SsSk_pooled_tajimad/"
echo ""
echo "Key output files (prefix = run tag):"
echo "  <tag>_pi.txt         per-window Pi per population"
echo "  <tag>_dxy.txt        per-window Dxy between population pairs"
echo "  <tag>_fst.txt        per-window Fst between population pairs"
echo "  <tag>_tajima_d.txt   per-window Tajima's D per population"
echo ""
echo "Summarize per-window values to per-population means for reporting."
