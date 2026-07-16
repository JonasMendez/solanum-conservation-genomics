#!/usr/bin/env Rscript
# =============================================================================
# 03_kinship_betadosage.R
# Estimates pairwise kinship (coancestry) between all sampled individuals
# using the beta-hat estimator of Weir & Goudet (2017), implemented via
# hierfstat::beta.dosage().
#
# Generalized to run for any of three taxa:
#   Ss — S. sandwicense (Oahu; F1 PAKA vs F2 PAHY)
#   Si — S. incompletum (Hawaii Island; wild vs nursery stock)
#   Sk — S. kavaiensis  (Kauai; F1 vs F2)
#
# ESTIMATOR: beta-hat (B^ij)
#   Measures kinship of each pair relative to the average kinship of ALL
#   pairs in the sample. Values < 0 = less related than average (prioritize
#   as cross candidates); values > 0 = more related than average.
#
# USAGE:
#   Rscript scripts/03_kinship_betadosage.R          # runs Ss (default)
#   Rscript scripts/03_kinship_betadosage.R Ss
#   Rscript scripts/03_kinship_betadosage.R Si
#   Rscript scripts/03_kinship_betadosage.R Sk
#
# INPUT VCFs (must be pre-filtered with 03_prefilter_vcf_popgen.sh):
#   Ss: A2_Ss_individual_noMAF_allSNP/populations.snps.filtered.vcf.gz
#   Si: A3_Si_individual_noMAF_allSNP/populations.snps.filtered.vcf.gz
#   Sk: A4_Sk_individual_noMAF_allSNP/populations.snps.filtered.vcf.gz
#
# NOTE: HWE filter is skipped in prefilter for individual-level popmaps
#   (each sample is its own population; HWE undefined at n=1).
#
# Requirements (conda env: popgen):
#   vcfR, hierfstat, dplyr, readr
#
# Run from project root.
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(hierfstat)
  library(dplyr)
  library(readr)
})

# =============================================================================
# ── DATASET SELECTION ─────────────────────────────────────────────────────────
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
TAXON <- if (length(args) >= 1) toupper(args[1]) else "Sk"

if (!TAXON %in% c("Ss", "Si", "Sk")) {
  stop("TAXON must be one of: Ss, Si, Sk\n",
       "Usage: Rscript scripts/03_kinship_betadosage.R [Ss|Si|Sk]")
}

message("Running kinship analysis for taxon: ", TAXON)

# =============================================================================
# ── PER-TAXON CONFIGURATION ───────────────────────────────────────────────────
# =============================================================================

POPRUN_BASE <- "02_datasets/populations_runs"
META_FILE   <- "metadata/Solanum_Metadata.csv"

# Configuration list for each taxon
configs <- list(

  Ss = list(
    taxon_label  = "Solanum sandwicense",
    vcf_tag      = "A2_Ss_individual_noMAF_allSNP",
    species_abbrev = "Ss",
    out_dir      = "03_analyses/kinship/Ss",
    prefix       = "Ss",
    ind_cr_min   = 0.75,   # lowered to retain Ss_PAKA_01 (callrate=0.799)
    # Metadata columns to annotate pairs with
    group_cols   = c("generation", "population_code"),
    # Function to build Cross_type label from group columns
    cross_type_fn = function(df) {
      dplyr::case_when(
        df$Gen1 == "F1" & df$Gen2 == "F1" ~ "F1xF1",
        df$Gen1 == "F2" & df$Gen2 == "F2" ~ "F2xF2",
        (df$Gen1 == "F1" & df$Gen2 == "F2") |
        (df$Gen1 == "F2" & df$Gen2 == "F1") ~ "F1xF2",
        TRUE ~ "unknown"
      )
    },
    # Metadata corrections specific to this taxon
    meta_corrections = function(meta) {
      meta %>% dplyr::mutate(
        generation      = dplyr::if_else(sample_id == "Ss_PAKA_10", "F1",      generation),
        population_code = dplyr::if_else(sample_id == "Ss_PAKA_10", "PAKA",    population_code),
        location        = dplyr::if_else(sample_id == "Ss_PAKA_10", "Palikea", location)
      )
    }
  ),

  Si = list(
    taxon_label  = "Solanum incompletum",
    vcf_tag      = "A3_Si_wildstock_noMAF_allSNP",
    species_abbrev = "Si",
    out_dir      = "03_analyses/kinship/Si",
    prefix       = "Si",
    ind_cr_min   = 0.80,
    group_cols   = c("generation", "location"),
    cross_type_fn = function(df) {
      dplyr::case_when(
        df$Gen1 == "wild"  & df$Gen2 == "wild"  ~ "wildxwild",
        df$Gen1 == "stock" & df$Gen2 == "stock" ~ "stockxstock",
        (df$Gen1 == "wild"  & df$Gen2 == "stock") |
        (df$Gen1 == "stock" & df$Gen2 == "wild")  ~ "wildxstock",
        TRUE ~ "unknown"
      )
    },
    meta_corrections = function(meta) {
      # Si_HAPAMA14_04 generation corrected from "unknown" to "stock"
      meta %>% dplyr::mutate(
        generation = dplyr::if_else(sample_id == "Si_HAPAMA14_04", "stock", generation)
      )
    }
  ),

  Sk = list(
    taxon_label  = "Solanum kavaiensis",
    vcf_tag      = "A4_Sk_F1F2_noMAF_allSNP",
    species_abbrev = "Sk",
    out_dir      = "03_analyses/kinship/Sk",
    prefix       = "Sk",
    ind_cr_min   = 0.80,
    group_cols   = c("generation", "population_code"),
    cross_type_fn = function(df) {
      dplyr::case_when(
        df$Gen1 == "F1" & df$Gen2 == "F1" ~ "F1xF1",
        df$Gen1 == "F2" & df$Gen2 == "F2" ~ "F2xF2",
        (df$Gen1 == "F1" & df$Gen2 == "F2") |
        (df$Gen1 == "F2" & df$Gen2 == "F1") ~ "F1xF2",
        TRUE ~ "unknown"
      )
    },
    meta_corrections = function(meta) meta   # no corrections needed for Sk
  )

)

cfg <- configs[[TAXON]]

# Global thresholds (same for all taxa)
LOCUS_CR_MIN          <- 0.80
IND_CR_MIN            <- cfg$ind_cr_min
KINSHIP_FLAG_THRESHOLD <- 0.0
TOP_N_CANDIDATES      <- 10

VCF_PATH <- file.path(POPRUN_BASE, cfg$vcf_tag, "populations.snps.filtered.vcf.gz")
OUT_DIR  <- cfg$out_dir
PREFIX   <- cfg$prefix

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================

is_missing_gt <- function(gt) {
  is.na(gt) | gt %in% c("./.", ".|.", ".", "..")
}

gt_to_dosage <- function(gt_matrix) {
  dos <- matrix(NA_real_, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix),
                dimnames = dimnames(gt_matrix))
  dos[gt_matrix == "0/0" | gt_matrix == "0|0"] <- 0
  dos[gt_matrix == "0/1" | gt_matrix == "1/0" |
      gt_matrix == "0|1" | gt_matrix == "1|0"] <- 1
  dos[gt_matrix == "1/1" | gt_matrix == "1|1"] <- 2
  dos
}

# =============================================================================
# ── MAIN ──────────────────────────────────────────────────────────────────────
# =============================================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load and filter metadata ──────────────────────────────────────────────────
message("Reading metadata: ", META_FILE)
meta <- readr::read_csv(META_FILE, show_col_types = FALSE) %>%
  dplyr::filter(species_abbrev == cfg$species_abbrev) %>%
  dplyr::select(sample_id, species_abbrev, generation, population_code, location)

# Apply taxon-specific metadata corrections
meta <- cfg$meta_corrections(meta)
message("  Samples in metadata: ", nrow(meta))

# ── Load VCF ──────────────────────────────────────────────────────────────────
message("\nReading VCF: ", VCF_PATH)
if (!file.exists(VCF_PATH)) {
  stop("VCF not found: ", VCF_PATH,
       "\nRun 03_prefilter_vcf_popgen.sh first to generate the filtered VCF.")
}

vcf <- read.vcfR(VCF_PATH, verbose = FALSE)
message("  Raw SNPs : ", nrow(vcf))
message("  Samples  : ", ncol(vcf@gt) - 1)

# ── Extract GT matrix ─────────────────────────────────────────────────────────
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# ── Locus call-rate filter ────────────────────────────────────────────────────
locus_cr  <- apply(gt, 1, function(x) mean(!is_missing_gt(x)))
keep_loci <- locus_cr >= LOCUS_CR_MIN
gt        <- gt[keep_loci, , drop = FALSE]
message("  SNPs after locus call-rate filter (>=", LOCUS_CR_MIN, "): ", nrow(gt))

# ── Individual call-rate check ────────────────────────────────────────────────
ind_cr <- apply(gt, 2, function(x) mean(!is_missing_gt(x)))

cr_df <- data.frame(
  sample_id     = names(ind_cr),
  ind_callrate  = as.numeric(ind_cr),
  passed_filter = ind_cr >= IND_CR_MIN
) %>%
  dplyr::left_join(meta, by = "sample_id") %>%
  dplyr::arrange(ind_callrate)

write.csv(cr_df,
          file.path(OUT_DIR, paste0(PREFIX, "_individual_callrates.csv")),
          row.names = FALSE)

low_cr <- cr_df %>% dplyr::filter(!passed_filter)
if (nrow(low_cr) > 0) {
  message("  WARNING: ", nrow(low_cr), " individual(s) below call-rate threshold (",
          IND_CR_MIN, "):")
  for (i in seq_len(nrow(low_cr))) {
    message("    ", low_cr$sample_id[i],
            "  callrate=", round(low_cr$ind_callrate[i], 4))
  }
  message("  Retained but flagged. Raise IND_CR_MIN to exclude.")
} else {
  message("  All individuals pass call-rate threshold (>=", IND_CR_MIN, ")")
}

keep_inds <- ind_cr >= IND_CR_MIN
gt        <- gt[, keep_inds, drop = FALSE]
message("  Samples retained for kinship: ", ncol(gt))

# ── Convert GT to dosage and transpose ────────────────────────────────────────
message("\nConverting genotypes to allele dosage (0/1/2)...")
dos <- gt_to_dosage(gt)

# beta.dosage expects individuals x loci; extract.gt returns loci x individuals
dos <- t(dos)
message("  Dosage matrix: ", nrow(dos), " individuals x ", ncol(dos), " loci")

# ── Compute pairwise kinship ──────────────────────────────────────────────────
message("\nComputing pairwise kinship (hierfstat::beta.dosage)...")

beta_mat <- beta.dosage(dos, inb = FALSE)

message("  Individuals passed in  : ", nrow(dos))
message("  Beta matrix dimensions : ", nrow(beta_mat), " x ", ncol(beta_mat))

# Handle potential internal sample dropping by beta.dosage
if (nrow(beta_mat) != nrow(dos)) {
  message("  NOTE: beta.dosage dropped ", nrow(dos) - nrow(beta_mat),
          " individual(s) with insufficient data.")
  if (!is.null(rownames(beta_mat))) {
    retained_samples <- rownames(beta_mat)
    message("  Retained: ", paste(retained_samples, collapse = ", "))
  } else {
    retained_samples <- rownames(dos)[seq_len(nrow(beta_mat))]
    message("  WARNING: beta.dosage returned no rownames — assuming first ",
            nrow(beta_mat), " individuals retained. Verify results carefully.")
  }
} else {
  retained_samples <- rownames(dos)
}

rownames(beta_mat) <- retained_samples
colnames(beta_mat) <- retained_samples

message("  Kinship matrix: ", nrow(beta_mat), " x ", ncol(beta_mat))
message("  Beta range    : [", round(min(beta_mat, na.rm = TRUE), 4),
        ", ", round(max(beta_mat, na.rm = TRUE), 4), "]")

# ── Save full pairwise matrix ─────────────────────────────────────────────────
mat_df <- as.data.frame(beta_mat)
mat_df <- cbind(sample_id = rownames(mat_df), mat_df)
write.csv(mat_df,
          file.path(OUT_DIR, paste0(PREFIX, "_kinship_beta_matrix.csv")),
          row.names = FALSE)
message("\nSaved: ", PREFIX, "_kinship_beta_matrix.csv")

# ── Build ranked pairs table ──────────────────────────────────────────────────
message("Building ranked pairs table...")

samps <- rownames(beta_mat)
pairs <- do.call(rbind, lapply(seq_len(nrow(beta_mat) - 1), function(i) {
  do.call(rbind, lapply(seq(i + 1, ncol(beta_mat)), function(j) {
    data.frame(Ind1 = samps[i], Ind2 = samps[j],
               Beta_ij = beta_mat[i, j], stringsAsFactors = FALSE)
  }))
}))

# Annotate pairs with metadata (using first group column as "Gen" for cross type)
meta_slim <- meta %>%
  dplyr::select(sample_id, generation, population_code)

pairs <- pairs %>%
  dplyr::left_join(meta_slim, by = c("Ind1" = "sample_id")) %>%
  dplyr::rename(Gen1 = generation, Pop1 = population_code) %>%
  dplyr::left_join(meta_slim, by = c("Ind2" = "sample_id")) %>%
  dplyr::rename(Gen2 = generation, Pop2 = population_code) %>%
  dplyr::mutate(
    Cross_type      = cfg$cross_type_fn(.),
    Cross_candidate = Beta_ij < KINSHIP_FLAG_THRESHOLD
  ) %>%
  dplyr::arrange(Beta_ij)

write.csv(pairs,
          file.path(OUT_DIR, paste0(PREFIX, "_kinship_pairs_ranked.csv")),
          row.names = FALSE)
message("Saved: ", PREFIX, "_kinship_pairs_ranked.csv  (", nrow(pairs), " pairs)")

# ── Save cross candidates ─────────────────────────────────────────────────────
candidates <- pairs %>% dplyr::filter(Cross_candidate)
write.csv(candidates,
          file.path(OUT_DIR, paste0(PREFIX, "_kinship_cross_candidates.csv")),
          row.names = FALSE)
message("Saved: ", PREFIX, "_kinship_cross_candidates.csv  (",
        nrow(candidates), " candidate pairs with Beta < ", KINSHIP_FLAG_THRESHOLD, ")")

# ── Per-individual mean kinship ───────────────────────────────────────────────
diag(beta_mat) <- NA
ind_mean_beta <- rowMeans(beta_mat, na.rm = TRUE)

ind_summary <- data.frame(
  sample_id           = names(ind_mean_beta),
  mean_beta_to_others = as.numeric(ind_mean_beta)
) %>%
  dplyr::left_join(meta, by = "sample_id") %>%
  dplyr::arrange(mean_beta_to_others)

write.csv(ind_summary,
          file.path(OUT_DIR, paste0(PREFIX, "_kinship_individual_mean.csv")),
          row.names = FALSE)
message("Saved: ", PREFIX, "_kinship_individual_mean.csv")

# ── Console summary ───────────────────────────────────────────────────────────
message("\n", strrep("=", 62))
message("KINSHIP SUMMARY — ", cfg$taxon_label)
message(strrep("=", 62))
message(sprintf("  Taxon              : %s (%s)", cfg$taxon_label, TAXON))
message(sprintf("  SNPs used          : %d", ncol(dos)))
message(sprintf("  Individuals        : %d passed in, %d retained in matrix",
                nrow(dos), nrow(beta_mat)))
message(sprintf("  Total pairs        : %d", nrow(pairs)))
message(sprintf("  Cross candidates   : %d  (Beta < %.2f)",
                nrow(candidates), KINSHIP_FLAG_THRESHOLD))
message(sprintf("  Beta range         : [%.4f, %.4f]",
                min(pairs$Beta_ij, na.rm = TRUE),
                max(pairs$Beta_ij, na.rm = TRUE)))
message(sprintf("  Mean pairwise Beta : %.4f", mean(pairs$Beta_ij, na.rm = TRUE)))

message("\nTop ", TOP_N_CANDIDATES, " most genetically distinct pairs (recommended crosses):")
top <- head(pairs, TOP_N_CANDIDATES)
for (i in seq_len(nrow(top))) {
  message(sprintf("  %2d. %-20s x %-20s  Beta=%.4f  [%s]",
                  i, top$Ind1[i], top$Ind2[i],
                  top$Beta_ij[i], top$Cross_type[i]))
}

message("\nMost genetically unique individuals (lowest mean kinship to others):")
top_ind <- head(ind_summary, 5)
for (i in seq_len(nrow(top_ind))) {
  message(sprintf("  %2d. %-20s  mean_Beta=%.4f  [%s / %s]",
                  i, top_ind$sample_id[i],
                  top_ind$mean_beta_to_others[i],
                  top_ind$generation[i],
                  top_ind$population_code[i]))
}

message("\n", strrep("=", 62))
message("All outputs in: ", OUT_DIR)
message("\nINTERPRETATION GUIDE:")
message("  Beta < 0  : pair less related than sample average -> PRIORITIZE for crosses")
message("  Beta ~ 0  : pair has average relatedness")
message("  Beta > 0  : pair more related than average -> DEPRIORITIZE")
message("  Beta ~0.25: ~first cousins (approximate, relative to sample average)")
message("  Beta ~0.5 : ~full siblings or parent-offspring (approximate)")
message("")
message("  Individual mean_beta: lowest values = most genetically unique.")
message("  Prioritize these individuals for inclusion in any cross.")
