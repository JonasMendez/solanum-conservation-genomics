#!/usr/bin/env bash
# =============================================================================
# 02_populations_runs.sh
# Runs all required `populations` exports from the single gstacks assembly.
#
# Each run produces a distinct dataset for a specific analysis question.
# Outputs go to 02_datasets/populations_runs/<tag>/
#
# Stacks locally installed: ~/local/stacks-2.68-noavx/bin
# Run from project root.
#
# Summary of runs and their purposes:
# ---------------------------------------------------------------------------
# TAG                                   POPMAP                  PURPOSE
# ---------------------------------------------------------------------------
# A1_all_taxa_MAF_randomSNP             all_taxa_species        PCA, ADMIXTURE (all taxa)
# A1_all_taxa_noMAF_randomSNP           all_taxa_species        He/Ho/Dxy/Fst via pixy (VCF only)
# A1_all_taxa_noMAF_allSNP              all_taxa_species        Pi/TajD/private alleles; GADMA2 input
# A1_all_taxa_phylo                     all_taxa_individual     IQ-TREE (--fasta-loci + --phylip-var-all)
# A1_all_taxa_indiv_phylip              all_taxa_individual     Per-individual PHYLIP → IUPAC VCF for pixy
# A2_Ss_F1F2_noMAF_randomSNP           Ss_F1F2                 PCA/ADMIXTURE/He/Ho/FIS (VCF only)
# A2_Ss_F1F2_noMAF_allSNP              Ss_F1F2                 Pi/private alleles (Ss)
# A2_Ss_individual_noMAF_allSNP        Ss_individual           Per-individual private alleles/kinship
# A2_Ss_indiv_phylip                   Ss_individual           Per-individual PHYLIP → IUPAC VCF for pixy
# A3_Si_wildstock_noMAF_randomSNP      Si_wildstock            PCA/ADMIXTURE/He/Ho/Dxy/Fst (VCF only)
# A3_Si_wildstock_noMAF_allSNP         Si_wildstock            Pi/TajD/private alleles (Si)
# A3_Si_individual_noMAF_allSNP        Si_individual           Per-individual private alleles/kinship
# A3_Si_indiv_phylip                   Si_individual           Per-individual PHYLIP → IUPAC VCF for pixy
# A4_Sk_F1F2_noMAF_randomSNP           Sk_F1F2                 PCA/ADMIXTURE/He/Ho/Dxy/Fst (VCF only)
# A4_Sk_F1F2_noMAF_allSNP              Sk_F1F2                 Pi/TajD/private alleles (Sk)
# A4_Sk_individual_noMAF_allSNP        Sk_individual           Per-individual private alleles/kinship
# A4_Sk_indiv_phylip                   Sk_individual           Per-individual PHYLIP → IUPAC VCF for pixy
# ---------------------------------------------------------------------------
# PHYLIP notes:
#   --phylip-var-all on GROUP-level popmaps produces one consensus sequence per
#   population group — NOT suitable for IUPAC VCF conversion (pixy needs one
#   sequence per individual). Therefore the four pixy-target datasets each have
#   a dedicated individual-level PHYLIP run (_indiv_phylip tags above).
#   These PHYLIP outputs are post-processed by 02_fix_phylip.sh (de-interleaving
#   + label fixing) before being passed to the IUPAC VCF converter.
#
#   --fasta-loci is used ONLY for A1_all_taxa_phylo (IQ-TREE supermatrix
#   concatenation via 03b_concat_stacks_fasta_loci_to_phylip.py).
#
# NOTE: A1_all_taxa_noMAF_allSNP VCF is also used as GADMA2/easySFS input (A5).
# =============================================================================
set -euo pipefail

export PATH="$HOME/local/stacks-2.68-noavx/bin:$PATH"

# ── global settings ────────────────────────────────────────────────────────
RUN_DIR="01_stacks/runs/run_m3_M2_n1"
OUTBASE="02_datasets/populations_runs"
LOGDIR="02_datasets/logs"
THREADS=8

R_GLOBAL=0.70   # between-taxa runs (-R): 70% of ALL samples must be genotyped
R_TAXON=0.70    # within-taxon runs (-r): 70% within each population group

mkdir -p "${OUTBASE}" "${LOGDIR}"

[[ -d "${RUN_DIR}" ]] || { echo "ERROR: Stacks run not found: ${RUN_DIR}"; exit 1; }

# ── helper ────────────────────────────────────────────────────────────────
run_populations() {
    local tag="$1"
    local popmap="$2"
    shift 2
    local extra_args=("$@")

    local outdir="${OUTBASE}/${tag}"
    local logfile="${LOGDIR}/populations_${tag}.log"
    mkdir -p "${outdir}"

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  TAG    : ${tag}"
    echo "  POPMAP : ${popmap}"
    echo "  ARGS   : ${extra_args[*]}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    populations \
        -P "${RUN_DIR}" \
        -M "${popmap}" \
        -O "${outdir}" \
        -t "${THREADS}" \
        "${extra_args[@]}" \
        2>&1 | tee "${logfile}"

    echo "  Done → ${outdir}"
    echo ""
}

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — All three taxa (Ss + Sk + Si)
# Global -R: 70% of all samples must be genotyped at a locus
# ══════════════════════════════════════════════════════════════════════════════

# PCA + ADMIXTURE: MAF filtered, one random SNP per locus
run_populations \
    "A1_all_taxa_MAF_randomSNP" \
    "metadata/popmaps/popmap_all_taxa_species.tsv" \
    -R "${R_GLOBAL}" \
    --min-maf 0.05 \
    --write-random-snp \
    --vcf

# pixy (He/Ho/Dxy/Fst): no MAF, one random SNP per locus, VCF only
# PHYLIP for pixy is generated separately by A1_all_taxa_indiv_phylip below
run_populations \
    "A1_all_taxa_noMAF_randomSNP" \
    "metadata/popmaps/popmap_all_taxa_species.tsv" \
    -R "${R_GLOBAL}" \
    --write-random-snp \
    --vcf

# Pi / Tajima's D / private alleles + GADMA2/easySFS input: no MAF, all SNPs
run_populations \
    "A1_all_taxa_noMAF_allSNP" \
    "metadata/popmaps/popmap_all_taxa_species.tsv" \
    -R "${R_GLOBAL}" \
    --vcf

# IQ-TREE supermatrix: individual-level popmap, per-locus FASTA + PHYLIP
run_populations \
    "A1_all_taxa_phylo" \
    "metadata/popmaps/popmap_all_taxa_individual.tsv" \
    -R "${R_GLOBAL}" \
    --fasta-loci \
    --phylip-var-all

# pixy PHYLIP source: individual-level popmap, PHYLIP only (no --fasta-loci)
# Produces one sequence per individual → correct input for IUPAC VCF converter
run_populations \
    "A1_all_taxa_indiv_phylip" \
    "metadata/popmaps/popmap_all_taxa_individual.tsv" \
    -R "${R_GLOBAL}" \
    --phylip-var-all

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — Ss: Palikea F1 (PAKA) vs Pahole F2 (PAHY)
# Within-taxon: -p 1 -r 0.70
# ══════════════════════════════════════════════════════════════════════════════

# PCA + ADMIXTURE + pixy (He/Ho/FIS): no MAF, one random SNP per locus, VCF only
run_populations \
    "A2_Ss_F1F2_noMAF_randomSNP" \
    "metadata/popmaps/popmap_Ss_F1F2.tsv" \
    -p 1 -r "${R_TAXON}" \
    --write-random-snp \
    --vcf

# Pi + private alleles: no MAF, all SNPs per locus
run_populations \
    "A2_Ss_F1F2_noMAF_allSNP" \
    "metadata/popmaps/popmap_Ss_F1F2.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# Per-individual private alleles + kinship input
run_populations \
    "A2_Ss_individual_noMAF_allSNP" \
    "metadata/popmaps/popmap_Ss_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# pixy PHYLIP source: individual-level popmap, PHYLIP only
run_populations \
    "A2_Ss_indiv_phylip" \
    "metadata/popmaps/popmap_Ss_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --phylip-var-all

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — Si: wild vs nursery stock
# Within-taxon: -p 1 -r 0.70
# ══════════════════════════════════════════════════════════════════════════════

# PCA + ADMIXTURE + pixy (He/Ho/Dxy/Fst): no MAF, one random SNP, VCF only
run_populations \
    "A3_Si_wildstock_noMAF_randomSNP" \
    "metadata/popmaps/popmap_Si_wildstock.tsv" \
    -p 1 -r "${R_TAXON}" \
    --write-random-snp \
    --vcf

# Pi + Tajima's D + private alleles: no MAF, all SNPs per locus
run_populations \
    "A3_Si_wildstock_noMAF_allSNP" \
    "metadata/popmaps/popmap_Si_wildstock.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# Per-individual private alleles + kinship input
run_populations \
    "A3_Si_individual_noMAF_allSNP" \
    "metadata/popmaps/popmap_Si_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# pixy PHYLIP source: individual-level popmap, PHYLIP only
run_populations \
    "A3_Si_indiv_phylip" \
    "metadata/popmaps/popmap_Si_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --phylip-var-all

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4 — Sk: F1 vs F2
# Within-taxon: -p 1 -r 0.70
# Note: Sk F1 has only 3 samples — keep an eye on SNP recovery
# ══════════════════════════════════════════════════════════════════════════════

# PCA + ADMIXTURE + pixy (He/Ho/Dxy/Fst): no MAF, one random SNP, VCF only
run_populations \
    "A4_Sk_F1F2_noMAF_randomSNP" \
    "metadata/popmaps/popmap_Sk_F1F2.tsv" \
    -p 1 -r "${R_TAXON}" \
    --write-random-snp \
    --vcf

# Pi + Tajima's D + private alleles: no MAF, all SNPs per locus
run_populations \
    "A4_Sk_F1F2_noMAF_allSNP" \
    "metadata/popmaps/popmap_Sk_F1F2.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# Per-individual private alleles + kinship input
run_populations \
    "A4_Sk_individual_noMAF_allSNP" \
    "metadata/popmaps/popmap_Sk_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --vcf

# pixy PHYLIP source: individual-level popmap, PHYLIP only
run_populations \
    "A4_Sk_indiv_phylip" \
    "metadata/popmaps/popmap_Sk_individual.tsv" \
    -p 1 -r "${R_TAXON}" \
    --phylip-var-all

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "All populations runs complete."
echo "Outputs in: ${OUTBASE}/"
echo ""
echo "PHYLIP files for de-interleaving + label fix (02_fix_phylip.sh):"
echo "  ${OUTBASE}/A1_all_taxa_indiv_phylip/populations.all.phylip"
echo "  ${OUTBASE}/A2_Ss_indiv_phylip/populations.all.phylip"
echo "  ${OUTBASE}/A3_Si_indiv_phylip/populations.all.phylip"
echo "  ${OUTBASE}/A4_Sk_indiv_phylip/populations.all.phylip"
echo ""
echo "IQ-TREE inputs (after 02_fix_phylip.sh + supermatrix concat):"
echo "  ${OUTBASE}/A1_all_taxa_phylo/  (--fasta-loci + --phylip-var-all)"
echo ""
echo "Next steps:"
echo "  bash scripts/02_prep_vcfs.sh"
echo "  bash scripts/02_fix_phylip.sh"
echo "  bash scripts/02_missingness_report.sh"
