#!/usr/bin/env bash
# =============================================================================
# 04_easysfs.sh
# Generates the joint site frequency spectrum (jSFS) for GADMA2 demographic
# inference using easySFS.
#
# WORKFLOW:
#   Step 1 (PREVIEW): Run easySFS --preview to assess the tradeoff between
#     projection size and number of segregating sites for each population.
#     Review the preview output to choose projection values before Step 2.
#
#   Step 2 (GENERATE): Run easySFS with chosen projections to produce the
#     jSFS input files for GADMA2.
#
# INPUT:
#   VCF : A1_all_taxa_noMAF_allSNP paralog-filtered VCF
#         (populations.snps.filtered.vcf.gz)
#   POPMAP : popmap_all_taxa_species.tsv (Ss/Sk/Si species-level groups)
#
# KEY FLAGS:
#   -a  : keep ALL SNPs per RAD locus (not one random SNP). This maximizes
#         segregating sites for the SFS. Linkage within 36bp RAD loci is
#         negligible for SFS-based demographic inference.
#   --order Si,Ss,Sk : population order matching GADMA2 model topology
#         ((Ss,Sk),Si) with outgroup Si listed first. GADMA2 and dadi/moments
#         expect consistent ordering throughout.
#   --total-length : estimated total sequence length = n_loci x 36bp.
#         Needed for accurate monomorphic (zero-bin) site count in the SFS.
#         Extracted from the Stacks populations log below.
#   Folded SFS: --unfolded flag is NOT used. Without outgroup polarization
#         data, the folded (minor allele frequency) SFS is appropriate.
#
# PROJECTION CHOICE GUIDANCE (review Step 1 output):
#   For each population, plot or inspect the segregating sites vs. projection
#   size table. Choose the projection at the "elbow" where additional
#   haplotypes yield diminishing returns in segregating sites. This is
#   typically where the curve flattens.
#   - Maximum possible projection = 2 x n_samples (diploid haplotypes)
#     Ss: max = 2x15 = 30  Sk: max = 2x8 = 16  Si: max = 2x14 = 28
#   - Projections do NOT need to be equal across populations.
#   - Smaller projections = fewer sites but less missing data bias.
#   - After reviewing preview, set PROJ_SS, PROJ_SK, PROJ_SI below and
#     run Step 2.
#
# Requirements (conda env: easysfs or gadma):
#   easySFS.py
#
# Run from project root.
# =============================================================================
set -euo pipefail

VCF="02_datasets/populations_runs/A1_all_taxa_noMAF_allSNP/populations.snps.filtered.vcf.gz"
POPMAP="metadata/popmaps/popmap_all_taxa_species.tsv"
OUTBASE="03_analyses/gadma"
PREVIEW_DIR="${OUTBASE}/easysfs_preview"
SFS_DIR="${OUTBASE}/easysfs_sfs"
LOGDIR="${OUTBASE}/logs"

mkdir -p "${PREVIEW_DIR}" "${SFS_DIR}" "${LOGDIR}"

# ── Total sequence length ─────────────────────────────────────────────────
# Derived from the all-taxa supermatrix PHYLIP alignment (A1_all_taxa_phylo),
# which reported a total concatenated alignment length of 462,971 bp across
# all retained RAD loci. This value is used as the --total-length argument
# for easySFS to accurately populate the monomorphic (zero-bin) site count.
TOTAL_LENGTH=462971
echo "Total sequence length: ${TOTAL_LENGTH} bp  (from supermatrix PHYLIP alignment)"
echo ""

# =============================================================================
# STEP 1 — PREVIEW: assess projection size vs segregating sites tradeoff
# =============================================================================
# echo "================================================================"
# echo "STEP 1: easySFS preview"
# echo "================================================================"
# echo ""
# echo "Preview output shows segregating sites for each possible projection"
# echo "value per population. Review this output to choose PROJ_SS, PROJ_SK,"
# echo "PROJ_SI before running Step 2."
# echo ""
# echo "Projection selection guide:"
# echo "  Max haplotypes: Ss=30 (n=15), Sk=16 (n=8), Si=28 (n=14)"
# echo "  Choose the projection at the elbow of the sites-vs-haplotypes curve."
# echo "  Projections need not be equal across populations."
# echo ""

# easySFS.py \
#     -i "${VCF}" \
#     -p "${POPMAP}" \
#     --preview \
#     --order Si,Ss,Sk \
#     --ploidy 2 \
#     -a \
#     -y \
#     2>&1 | tee "${PREVIEW_DIR}/easysfs_preview.txt"

# echo ""
# echo "Preview output saved: ${PREVIEW_DIR}/easysfs_preview.txt"
# echo ""
# echo "================================================================"
# echo "REVIEW PREVIEW OUTPUT ABOVE THEN:"
# echo "  1. Set PROJ_SI, PROJ_SS, PROJ_SK in Step 2 below"
# echo "  2. Uncomment and run the Step 2 block"
# echo "================================================================"

# =============================================================================
# STEP 2 — GENERATE SFS: run after reviewing preview and setting projections
# =============================================================================
# Uncomment this block after reviewing the Step 1 preview output.
# Set projection values based on the elbow of the segregating sites curve.
# Order must match --order Si,Ss,Sk

PROJ_SI=16   # Set after reviewing preview (max=28)
PROJ_SS=30   # Set after reviewing preview (max=30)
PROJ_SK=12   # Set after reviewing preview (max=16)

# echo "================================================================"
# echo "STEP 2: Generating jSFS with projections Si=${PROJ_SI}, Ss=${PROJ_SS}, Sk=${PROJ_SK}"
# echo "================================================================"
# echo ""

# easySFS.py \
#     -i "${VCF}" \
#     -p "${POPMAP}" \
#     -o "${SFS_DIR}" \
#     --proj "${PROJ_SI},${PROJ_SS},${PROJ_SK}" \
#     --order Si,Ss,Sk \
#     --prefix "Solanum_3pop" \
#     --ploidy 2 \
#     --total-length "${TOTAL_LENGTH}" \
#     --dtype float \
#     -a \
#     -f \
#     -y \
#     2>&1 | tee "${LOGDIR}/easysfs_generate.log"

# echo ""
# echo "SFS files written to: ${SFS_DIR}/"
# echo ""
# echo "Key output files:"
# echo "  ${SFS_DIR}/dadi/Solanum_3pop-DSFS.fs  (dadi/moments format -- used by GADMA2)"
# echo ""
# echo "Next step: bash scripts/04_gadma2_run.sh"

#Step 3 — Testing alternative topoligy with Sk outgroup; filled in the variables manually:
# SFS_DIR_B="03_analyses/gadma/easysfs_sfs_Sk_outgroup"
# mkdir -p "${SFS_DIR_B}"

# easySFS.py \
#     -i 02_datasets/populations_runs/A1_all_taxa_noMAF_allSNP/populations.snps.filtered.vcf.gz \
#     -p metadata/popmaps/popmap_all_taxa_species.tsv \
#     -o "${SFS_DIR_B}" \
#     --proj 12,16,30 \
#     --order Sk,Si,Ss \
#     --prefix "Solanum_3pop_Sk_outgroup" \
#     --ploidy 2 \
#     --total-length 462971 \
#     --dtype float \
#     -a \
#     -f \
#     -y

