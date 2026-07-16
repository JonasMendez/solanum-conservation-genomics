#!/usr/bin/env Rscript
# =============================================================================
# 03_popgen_diversity_fst.R
# Population genetic diversity and differentiation for Hawaiian Solanum 2bRAD data.
#
# Computes per-dataset:
#   (A) Within-group diversity : Ho, He, FIS (hierfstat::basic.stats)
#   (B) Allelic richness       : rarefied (hierfstat::allelic.richness)
#   (C) Pairwise FST           : Weir & Cockerham (hierfstat::pairwise.WCfst)
#   (D) AMOVA                  : among groups (poppr::poppr.amova + pegas::randtest)
#
# FILTERING RATIONALE:
#   All analyses run from the no-MAF, all-SNP VCFs. MAF filtering is NOT applied
#   for He/Ho/FIS (rare alleles are genuine biological signal that contribute to
#   heterozygosity estimates), and is also not applied for FST/AMOVA because
#   rare private alleles are precisely the signal driving inter-taxon/group
#   differentiation — removing them would underestimate true FST. A call-rate
#   filter (locus genotyped in ≥80% of samples; individual genotyped at ≥80% of
#   loci) is applied to ensure robust genotype data without MAF-based censoring.
#   This is internally consistent: all statistics are computed on the same SNP set.
#
# DATASETS:
#   A1: All taxa (Ss + Sk + Si) — FST + AMOVA are the primary interest;
#       He/Ho/FIS computed per species as a bonus
#   A2: Ss F1 (PAKA) vs F2 (PAHY) — He/Ho/FIS primary; FST/AMOVA included
#   A3: Si wild vs stock           — He/Ho/FIS primary; FST/AMOVA included
#   A4: Sk F1 vs F2               — He/Ho/FIS primary; FST/AMOVA included
#
# VCF INPUTS: the all-SNP, no-MAF populations.snps.fixed.vcf.gz outputs
#   A1: 02_datasets/populations_runs/A1_all_taxa_noMAF_allSNP/
#   A2: 02_datasets/populations_runs/A2_Ss_F1F2_noMAF_allSNP/
#   A3: 02_datasets/populations_runs/A3_Si_wildstock_noMAF_allSNP/
#   A4: 02_datasets/populations_runs/A4_Sk_F1F2_noMAF_allSNP/
#
# POPMAP INPUTS: TSV popmaps from metadata/popmaps/
#
# Requirements (conda env: popgen or equivalent):
#   vcfR, adegenet, hierfstat, poppr, pegas, dplyr, tidyr, ggplot2, readr
#
# Usage:
#   Rscript scripts/03_popgen_diversity_fst.R
#
# Outputs written to: 03_analyses/popgen/<dataset_tag>/
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(poppr)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

set.seed(123)

# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================

POPRUN_BASE <- "02_datasets/populations_runs"
POPMAP_DIR  <- "metadata/popmaps"
OUT_BASE    <- "03_analyses/popgen"

# Call-rate filters applied to the all-SNP VCF within this script
LOCUS_CALLRATE_MIN <- 0.80   # locus genotyped in ≥80% of samples
IND_CALLRATE_MIN   <- 0.80   # individual genotyped at ≥80% of loci

# Allelic richness rarefaction: NULL = auto (min group n after filtering)
AR_MIN_N <- NULL

# AMOVA permutations
AMOVA_NREPET <- 999

# ── Dataset definitions ───────────────────────────────────────────────────
# Each entry: tag, vcf_tag (folder name), popmap file, group_col in popmap
DATASETS <- list(

  list(
    tag       = "A1_all_taxa",
    vcf_tag   = "A1_all_taxa_noMAF_allSNP",
    popmap    = "popmap_all_taxa_species.tsv",
    # Column 2 of the popmap = species_abbrev (Ss/Sk/Si)
    # FST + AMOVA are the primary outputs for this dataset
    primary   = "FST_AMOVA"
  ),

  list(
    tag       = "A2_Ss_F1F2",
    vcf_tag   = "A2_Ss_F1F2_noMAF_allSNP",
    popmap    = "popmap_Ss_F1F2.tsv",
    # Column 2 = population_code (PAKA/PAHY)
    primary   = "He_Ho_FIS"
  ),

  list(
    tag       = "A3_Si_wildstock",
    vcf_tag   = "A3_Si_wildstock_noMAF_allSNP",
    popmap    = "popmap_Si_wildstock.tsv",
    # Column 2 = generation (wild/stock)
    primary   = "He_Ho_FIS"
  ),

  list(
    tag       = "A4_Sk_F1F2",
    vcf_tag   = "A4_Sk_F1F2_noMAF_allSNP",
    popmap    = "popmap_Sk_F1F2.tsv",
    # Column 2 = generation (F1/F2)
    primary   = "He_Ho_FIS"
  )

)

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================

is_missing_gt <- function(gt) {
  is.na(gt) | gt %in% c("./.", ".|.", ".", "..")
}

read_popmap <- function(path) {
  if (!file.exists(path)) stop("Popmap not found: ", path)
  df <- readr::read_tsv(path, col_names = c("sample_id", "group"),
                         show_col_types = FALSE)
  df
}

load_and_filter_vcf <- function(vcf_path, popmap,
                                 locus_cr = LOCUS_CALLRATE_MIN,
                                 ind_cr   = IND_CALLRATE_MIN) {

  message("  Reading VCF: ", vcf_path)
  vcf <- read.vcfR(vcf_path, verbose = FALSE)

  # Subset to samples present in popmap (in VCF order)
  vcf_samples <- colnames(vcf@gt)[-1]
  keep_samp   <- vcf_samples %in% popmap$sample_id
  vcf <- vcf[, c(TRUE, keep_samp)]
  vcf_samples <- colnames(vcf@gt)[-1]
  message("  Samples matched to popmap: ", length(vcf_samples))

  # Extract GT matrix
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

  # Locus call-rate filter
  locus_cr_vec <- apply(gt, 1, function(x) mean(!is_missing_gt(x)))
  keep_loci    <- locus_cr_vec >= locus_cr
  vcf <- vcf[keep_loci, ]
  gt  <- gt[keep_loci, , drop = FALSE]
  message("  Loci after call-rate filter (>=", locus_cr, "): ", nrow(gt))

  # Individual call-rate filter
  ind_cr_vec <- apply(gt, 2, function(x) mean(!is_missing_gt(x)))
  keep_inds  <- ind_cr_vec >= ind_cr
  vcf  <- vcf[, c(TRUE, keep_inds)]
  gt   <- gt[, keep_inds, drop = FALSE]
  message("  Samples after call-rate filter (>=", ind_cr, "): ", ncol(gt))

  # Align popmap to filtered sample order
  pm_f <- popmap %>%
    filter(sample_id %in% colnames(gt)) %>%
    mutate(sample_id = factor(sample_id, levels = colnames(gt))) %>%
    arrange(sample_id) %>%
    mutate(sample_id = as.character(sample_id))

  list(vcf = vcf, gt = gt, popmap = pm_f,
       locus_cr = locus_cr_vec[keep_loci],
       ind_cr_all   = ind_cr_vec,
       ind_kept     = keep_inds)
}

build_genind <- function(vcf, popmap) {
  gi <- vcfR2genind(vcf)
  pop(gi) <- popmap$group
  adegenet::ploidy(gi) <- 2
  adegenet::strata(gi) <- data.frame(Group = pop(gi))
  gi <- adegenet::setPop(gi, ~Group)
  gi
}

# ── Analysis functions ────────────────────────────────────────────────────

run_diversity <- function(hfdat, genind, tag, out_dir, ar_min_n = NULL) {
  message("  Computing Ho / He / FIS (hierfstat::basic.stats)...")
  bs <- basic.stats(hfdat)

  Ho  <- colMeans(bs$Ho,  na.rm = TRUE)
  He  <- colMeans(bs$Hs,  na.rm = TRUE)
  Fis <- colMeans(bs$Fis, na.rm = TRUE)

  group_n <- as.numeric(table(pop(genind)))
  names(group_n) <- names(table(pop(genind)))

  div_df <- data.frame(
    Group    = names(Ho),
    N        = group_n[names(Ho)],
    Ho_mean  = as.numeric(Ho),
    He_mean  = as.numeric(He),
    Fis_mean = as.numeric(Fis)
  )

  # Allelic richness (rarefied to min group n)
  message("  Computing allelic richness (rarefied)...")
  min_n <- if (is.null(ar_min_n)) min(div_df$N) else ar_min_n
  ar    <- allelic.richness(hfdat, min.n = min_n)
  AR_mean <- colMeans(ar$Ar, na.rm = TRUE)
  div_df$AllelicRichness_mean <- as.numeric(AR_mean[div_df$Group])
  div_df$AR_rarefaction_n <- min_n

  out_path <- file.path(out_dir, paste0(tag, "_diversity_Ho_He_Fis_AR.csv"))
  write.csv(div_df, out_path, row.names = FALSE)
  message("  Saved: ", out_path)

  # Per-locus outputs for downstream plotting if needed
  ho_locus <- as.data.frame(bs$Ho)
  he_locus <- as.data.frame(bs$Hs)
  fis_locus <- as.data.frame(bs$Fis)
  write.csv(ho_locus,  file.path(out_dir, paste0(tag, "_Ho_perlocus.csv")))
  write.csv(he_locus,  file.path(out_dir, paste0(tag, "_He_perlocus.csv")))
  write.csv(fis_locus, file.path(out_dir, paste0(tag, "_Fis_perlocus.csv")))

  div_df
}

run_fst <- function(hfdat, tag, out_dir) {
  message("  Computing pairwise WC FST (hierfstat::pairwise.WCfst)...")
  fst_mat <- pairwise.WCfst(hfdat)

  # Matrix form
  out_mat <- file.path(out_dir, paste0(tag, "_pairwise_WCfst_matrix.csv"))
  write.csv(as.data.frame(fst_mat), out_mat, row.names = TRUE)
  message("  Saved: ", out_mat)

  # Long form (easier for plotting)
  fst_long <- as.data.frame(as.table(fst_mat)) %>%
    setNames(c("Group1", "Group2", "WC_FST")) %>%
    filter(!is.na(WC_FST), Group1 != Group2)
  out_long <- file.path(out_dir, paste0(tag, "_pairwise_WCfst_long.csv"))
  write.csv(fst_long, out_long, row.names = FALSE)

  fst_mat
}

run_amova <- function(genind, tag, out_dir, nrepet = AMOVA_NREPET) {
  message("  Running AMOVA (poppr::poppr.amova + ", nrepet, " permutations)...")

  # within=FALSE: avoids within-individual variance components that can be
  # unstable when genotype data has residual missing values
  amova_res  <- poppr.amova(genind, ~Group, within = FALSE)
  amova_test <- randtest(amova_res, nrepet = nrepet)

  # Variance components table
  comp <- amova_res$componentsofcovariance
  amova_tbl <- data.frame(
    Source  = rownames(comp),
    Sigma2  = comp[, 1],
    Percent = comp[, 1] / sum(comp[, 1]) * 100
  )
  out_comp <- file.path(out_dir, paste0(tag, "_AMOVA_components.csv"))
  write.csv(amova_tbl, out_comp, row.names = FALSE)
  message("  Saved: ", out_comp)

  # Permutation test result
  out_perm <- file.path(out_dir, paste0(tag, "_AMOVA_randtest.csv"))
  write.csv(
    data.frame(Observed = amova_test$obs, P_value = amova_test$pvalue),
    out_perm, row.names = FALSE
  )
  message("  AMOVA p-value: ", amova_test$pvalue)

  list(amova = amova_res, test = amova_test, components = amova_tbl)
}

# =============================================================================
# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
# =============================================================================

for (ds in DATASETS) {
  tag     <- ds$tag
  vcf_tag <- ds$vcf_tag
  primary <- ds$primary

  message("\n", strrep("━", 64))
  message("Dataset : ", tag, "  (primary: ", primary, ")")

  popmap_path <- file.path(POPMAP_DIR,  ds$popmap)
  out_dir     <- file.path(OUT_BASE, tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Select VCF: prefer HWE+het pre-filtered VCF if it exists.
  # Run 03_prefilter_vcf_popgen.sh first to generate the filtered VCF.
  # Falls back to fixed VCF with a warning if pre-filtering not yet done.
  vcf_filtered <- file.path(POPRUN_BASE, vcf_tag, "populations.snps.filtered.vcf.gz")
  vcf_fixed    <- file.path(POPRUN_BASE, vcf_tag, "populations.snps.fixed.vcf.gz")
  if (file.exists(vcf_filtered)) {
    vcf_path <- vcf_filtered
    message("  VCF: pre-filtered (HWE + max-het) -- ", basename(vcf_path))
  } else {
    vcf_path <- vcf_fixed
    message("  WARNING: pre-filtered VCF not found; using fixed VCF.")
    message("  Run 03_prefilter_vcf_popgen.sh first for reliable He/Ho/FIS results.")
  }

  # ── Validate inputs ───────────────────────────────────────────────────────
  if (!file.exists(vcf_path)) {
    message("SKIP — VCF not found: ", vcf_path); next
  }
  if (!file.exists(popmap_path)) {
    message("SKIP — popmap not found: ", popmap_path); next
  }

  # ── Load popmap ───────────────────────────────────────────────────────────
  popmap <- read_popmap(popmap_path)

  # ── Load and filter VCF ───────────────────────────────────────────────────
  dat <- load_and_filter_vcf(vcf_path, popmap)

  # Save per-individual call rates for QC reference
  # Report call rates for ALL samples (including those filtered out)
  # so excluded samples are visible with their actual call rate
  cr_df <- data.frame(
    sample_id    = names(dat$ind_cr_all),
    ind_callrate = as.numeric(dat$ind_cr_all),
    passed_filter = dat$ind_kept
  )
  # Add group from popmap where available
  cr_df <- merge(cr_df,
                 dat$popmap[, c("sample_id", "group")],
                 by = "sample_id", all.x = TRUE)
  cr_df <- cr_df[order(cr_df$sample_id), ]
  write.csv(cr_df, file.path(out_dir, paste0(tag, "_individual_callrates.csv")),
            row.names = FALSE)

  # ── Build genind + hierfstat objects ──────────────────────────────────────
  message("  Building genind object...")
  genind <- build_genind(dat$vcf, dat$popmap)
  hfdat  <- genind2hierfstat(genind)

  n_groups <- length(unique(dat$popmap$group))
  message("  Groups (n=", n_groups, "): ",
          paste(sort(unique(dat$popmap$group)), collapse = ", "))

  # ── (A) He / Ho / FIS / Allelic richness ─────────────────────────────────
  div_df <- run_diversity(hfdat, genind, tag, out_dir, ar_min_n = AR_MIN_N)
  message("  Diversity summary:")
  print(div_df, digits = 4)

  # ── (B) Pairwise FST ─────────────────────────────────────────────────────
  # Skip pairwise FST if only one group (shouldn't happen but safeguard)
  if (n_groups >= 2) {
    fst_mat <- run_fst(hfdat, tag, out_dir)
    message("  FST matrix:")
    print(round(fst_mat, 4))
  } else {
    message("  SKIP pairwise FST — only one group")
  }

  # ── (C) AMOVA ────────────────────────────────────────────────────────────
  if (n_groups >= 2) {
    amova_out <- run_amova(genind, tag, out_dir, nrepet = AMOVA_NREPET)
    message("  AMOVA variance components:")
    print(amova_out$components, digits = 4)
  } else {
    message("  SKIP AMOVA — only one group")
  }
}

message("\n", strrep("━", 64))
message("All popgen analyses complete.")
message("Outputs in: ", OUT_BASE)
message("")
message("Key output files per dataset:")
message("  <tag>_diversity_Ho_He_Fis_AR.csv   — mean He/Ho/FIS/AR per group")
message("  <tag>_pairwise_WCfst_matrix.csv    — WC FST matrix")
message("  <tag>_pairwise_WCfst_long.csv      — FST long format (for plotting)")
message("  <tag>_AMOVA_components.csv          — variance components + %")
message("  <tag>_AMOVA_randtest.csv            — permutation test p-value")
message("  <tag>_Ho_perlocus.csv / He / Fis   — per-locus values (for plotting distributions)")
message("  <tag>_individual_callrates.csv      — QC: per-sample genotyping rate")
