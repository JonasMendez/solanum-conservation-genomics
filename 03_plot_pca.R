#!/usr/bin/env Rscript
# =============================================================================
# 03_plot_pca.R
# Plots PCA results for all four Solanum datasets.
#
# Aesthetics are fully configurable per-dataset via the DATASET_CONFIG list:
#   - color_col   : metadata column to map to point color (e.g. "species_abbrev")
#   - shape_col   : metadata column to map to point shape (e.g. "generation")
#   - size_col    : metadata column to map to point SIZE (continuous, e.g. private alleles)
#                   OR NULL for constant size
#   - label_col   : metadata column to use as text label (NULL = no labels)
#   - pc_pairs    : list of PC axis pairs to plot, e.g. list(c("PC1","PC2"), c("PC1","PC3"))
#
# The script reads from metadata/Solanum_Metadata.csv by default but accepts
# an alternative metadata file path via the META_FILE variable below, as long
# as it contains a "sample_id" column plus whatever columns are referenced in
# DATASET_CONFIG. An optional SIZE_AUX_FILE can supply an additional per-sample
# numeric column (e.g. private allele counts) to use for point sizing; it must
# have columns "sample_id" and the column name referenced in size_col.
#
# Usage:
#   Rscript scripts/03_plot_pca.R
#   Rscript scripts/03_plot_pca.R --meta metadata/alt_metadata.csv
#   Rscript scripts/03_plot_pca.R --size_aux metadata/private_alleles_per_ind.tsv
#
# Outputs (per dataset, per PC pair):
#   04_plots/pca/<tag>/PCA_<tag>_<PCx>_<PCy>.png
#   04_plots/pca/<tag>/PCA_<tag>_<PCx>_<PCy>.pdf
#   04_plots/pca/<tag>/scree_<tag>.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# =============================================================================
# ── USER SETTINGS ─────────────────────────────────────────────────────────────
# =============================================================================

# Primary metadata file (must contain "sample_id" column)
META_FILE <- "metadata/Solanum_Metadata.csv"

# Optional auxiliary file with per-sample numeric values for point sizing
# (e.g. private allele counts). Must have "sample_id" + one numeric column.
# Set to NULL to disable; or override via --size_aux argument.
SIZE_AUX_FILE <- NULL

# PCA output base directory (one sub-folder per dataset tag)
PCA_BASE <- "03_analyses/pca"

# Plot output base directory
PLOT_BASE <- "04_plots/pca"

# Global aesthetics (overridable per-dataset in DATASET_CONFIG)
ALPHA           <- 0.90
POINT_SIZE_CONST <- 3.5    # used when no size_col is mapped
SIZE_RANGE      <- c(1, 10) # range when size_col is a continuous variable
LABEL_SIZE      <- 2.8
MAX_OVERLAPS    <- 80       # ggrepel max.overlaps

# ── Per-dataset configuration ──────────────────────────────────────────────
# Keys:
#   tag        : must match folder name under PCA_BASE/
#   color_col  : metadata column → point color  (NULL = single color)
#   shape_col  : metadata column → point shape  (NULL = single shape)
#   size_col   : metadata column → point size   (NULL = POINT_SIZE_CONST)
#                  if the column is numeric it scales continuously (SIZE_RANGE)
#                  if it is character it will be IGNORED with a warning
#   label_col  : metadata column → text label   (NULL = no labels)
#   pc_pairs   : list of character vectors, e.g. list(c("PC1","PC2"), c("PC1","PC3"))

DATASET_CONFIG <- list(

  A1_all_taxa_MAF_randomSNP = list(
    tag        = "A1_all_taxa_MAF_randomSNP",
    color_col  = "species_abbrev",   # Ss / Sk / Si
    shape_col  = "generation",       # F1 / F2 / wild / stock
    size_col   = NULL,               # constant size; swap to e.g. "private_alleles" if aux loaded
    label_col  = "sample_id",
    pc_pairs   = list(c("PC1","PC2"), c("PC1","PC3"), c("PC2","PC3"))
  ),

  A2_Ss_F1F2_noMAF_randomSNP = list(
    tag        = "A2_Ss_F1F2_noMAF_randomSNP",
    color_col  = "population_code",  # PAKA / PAHY
    shape_col  = "generation",       # F1 / F2
    size_col   = NULL,
    label_col  = "sample_id",
    pc_pairs   = list(c("PC1","PC2"), c("PC1","PC3"))
  ),

  A3_Si_wildstock_noMAF_randomSNP = list(
    tag        = "A3_Si_wildstock_noMAF_randomSNP",
    color_col  = "generation",       # wild / stock
    shape_col  = "location",         # site-level grouping
    size_col   = NULL,
    label_col  = "sample_id",
    pc_pairs   = list(c("PC1","PC2"), c("PC1","PC3"))
  ),

  A4_Sk_F1F2_noMAF_randomSNP = list(
    tag        = "A4_Sk_F1F2_noMAF_randomSNP",
    color_col  = "generation",       # F1 / F2
    shape_col  = "population_code",  # location-level grouping
    size_col   = NULL,
    label_col  = "sample_id",
    pc_pairs   = list(c("PC1","PC2"), c("PC1","PC3"))
  )

)

# =============================================================================
# ── COMMAND-LINE ARGUMENT PARSING ─────────────────────────────────────────────
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) && length(args) >= idx + 1) return(args[idx + 1])
  default
}

META_FILE     <- get_arg("--meta",     META_FILE)
SIZE_AUX_FILE <- get_arg("--size_aux", SIZE_AUX_FILE)

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================

read_meta <- function(path) {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  } else {
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
  }
}

read_size_aux <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  df <- if (ext == "csv") {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  } else {
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
  }
  if (!"sample_id" %in% colnames(df)) stop("SIZE_AUX_FILE must have a 'sample_id' column: ", path)
  df
}

read_plink_pca <- function(eigenvec_path, eigenval_path) {
  if (!file.exists(eigenvec_path)) stop("Missing eigenvec: ", eigenvec_path)
  if (!file.exists(eigenval_path)) stop("Missing eigenval: ", eigenval_path)

  ev <- suppressMessages(readr::read_table(
    eigenvec_path, col_names = FALSE, progress = FALSE, show_col_types = FALSE
  ))

  ncols <- ncol(ev)
  if (ncols < 4) stop("Unexpected eigenvec columns: ", ncols)

  pc_names <- paste0("PC", seq_len(ncols - 2))
  colnames(ev) <- c("FID", "IID", pc_names)

  evals <- suppressMessages(readr::read_table(
    eigenval_path, col_names = FALSE, progress = FALSE, show_col_types = FALSE
  )) |> dplyr::pull(1)

  pct <- evals / sum(evals) * 100

  list(scores = ev, eigenvals = evals, pct = pct, n_pcs = length(evals))
}

axis_label <- function(pc, pct_vec) {
  idx <- as.integer(str_remove(pc, "PC"))
  if (is.na(idx) || idx < 1 || idx > length(pct_vec)) return(pc)
  sprintf("%s (%.2f%%)", pc, pct_vec[idx])
}

save_plot <- function(p, path_prefix, width = 9, height = 7) {
  ggsave(paste0(path_prefix, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(path_prefix, ".pdf"), p, width = width, height = height)
  message("  Saved: ", basename(paste0(path_prefix, ".png")))
}

# =============================================================================
# ── PLOT FUNCTIONS ─────────────────────────────────────────────────────────────
# =============================================================================

plot_scree <- function(pca_obj, tag, out_dir) {
  n_show <- min(pca_obj$n_pcs, 15)
  df <- tibble(
    PC  = factor(paste0("PC", seq_len(n_show)), levels = paste0("PC", seq_len(n_show))),
    pct = pca_obj$pct[seq_len(n_show)]
  )

  p <- ggplot(df, aes(x = PC, y = pct)) +
    geom_col(fill = "#2E5496", width = 0.7) +
    geom_line(aes(group = 1), color = "#C00000", linewidth = 0.7) +
    geom_point(color = "#C00000", size = 2) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("Scree plot: ", tag),
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold")
    )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_plot(p, file.path(out_dir, paste0("scree_", tag)), width = 8, height = 5)
}

plot_pca_biplot <- function(pca_obj, cfg, meta, size_aux, pcx, pcy, out_dir) {

  df <- pca_obj$scores |>
    mutate(sample_id = IID) |>
    left_join(meta, by = "sample_id")

  # ── join size auxiliary data ─────────────────────────────────────────────
  use_size_col <- cfg$size_col
  if (!is.null(use_size_col) && !is.null(size_aux) && use_size_col %in% colnames(size_aux)) {
    df <- df |> left_join(size_aux |> select(sample_id, all_of(use_size_col)), by = "sample_id")
  } else if (!is.null(use_size_col) && !use_size_col %in% colnames(df)) {
    message("  NOTE: size_col '", use_size_col, "' not found in metadata or aux file — using constant size")
    use_size_col <- NULL
  }

  # Validate size_col is numeric if provided
  if (!is.null(use_size_col) && use_size_col %in% colnames(df)) {
    if (!is.numeric(df[[use_size_col]])) {
      message("  WARNING: size_col '", use_size_col, "' is not numeric — using constant size")
      use_size_col <- NULL
    }
  }

  # ── axis labels ──────────────────────────────────────────────────────────
  xlab <- axis_label(pcx, pca_obj$pct)
  ylab <- axis_label(pcy, pca_obj$pct)

  # ── base plot ────────────────────────────────────────────────────────────
  aes_base <- aes(x = .data[[pcx]], y = .data[[pcy]])

  p <- ggplot(df, aes_base) +
    theme_classic(base_size = 14) +
    labs(
      title  = paste0(cfg$tag, "  |  ", pcx, " vs ", pcy),
      x      = xlab,
      y      = ylab,
      color  = if (!is.null(cfg$color_col)) cfg$color_col else NULL,
      shape  = if (!is.null(cfg$shape_col)) cfg$shape_col else NULL,
      size   = if (!is.null(use_size_col))  use_size_col  else NULL
    ) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      legend.position = "right",
      legend.title    = element_text(size = 11),
      legend.text     = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4)

  # ── build aes for geom_point ─────────────────────────────────────────────
  point_aes <- list()
  if (!is.null(cfg$color_col) && cfg$color_col %in% colnames(df))
    point_aes$color <- rlang::sym(cfg$color_col)
  if (!is.null(cfg$shape_col) && cfg$shape_col %in% colnames(df))
    point_aes$shape <- rlang::sym(cfg$shape_col)

  if (!is.null(use_size_col) && use_size_col %in% colnames(df)) {
    # Samples with valid size values: sized by column
    df_has  <- df |> filter(!is.na(.data[[use_size_col]]))
    df_miss <- df |> filter(is.na(.data[[use_size_col]]))

    if (nrow(df_miss) > 0) {
      p <- p + do.call(
        geom_point,
        c(list(data = df_miss, mapping = do.call(aes, point_aes),
               size = POINT_SIZE_CONST, alpha = ALPHA))
      )
    }
    size_aes <- c(point_aes, list(size = rlang::sym(use_size_col)))
    p <- p +
      do.call(geom_point,
              c(list(data = df_has, mapping = do.call(aes, size_aes), alpha = ALPHA))) +
      scale_size_continuous(range = SIZE_RANGE)
  } else {
    p <- p + do.call(
      geom_point,
      c(list(mapping = do.call(aes, point_aes),
             size = POINT_SIZE_CONST, alpha = ALPHA))
    )
  }

  # ── labels ───────────────────────────────────────────────────────────────
  if (!is.null(cfg$label_col) && cfg$label_col %in% colnames(df)) {
    label_aes <- aes(label = .data[[cfg$label_col]])
    if (has_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        label_aes, size = LABEL_SIZE, max.overlaps = MAX_OVERLAPS,
        segment.color = "grey60", segment.size = 0.3
      )
    } else {
      p <- p + geom_text(label_aes, size = LABEL_SIZE, vjust = -0.7)
    }
  }

  # ── save ─────────────────────────────────────────────────────────────────
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(out_dir, paste0("PCA_", cfg$tag, "_", pcx, "_", pcy))
  save_plot(p, out_prefix)
}

# =============================================================================
# ── MAIN ──────────────────────────────────────────────────────────────────────
# =============================================================================

message("Reading metadata: ", META_FILE)
meta <- read_meta(META_FILE)
if (!"sample_id" %in% colnames(meta)) stop("Metadata must contain a 'sample_id' column")

size_aux <- NULL
if (!is.null(SIZE_AUX_FILE)) {
  message("Reading size auxiliary file: ", SIZE_AUX_FILE)
  size_aux <- read_size_aux(SIZE_AUX_FILE)
}

for (cfg in DATASET_CONFIG) {
  tag      <- cfg$tag
  pca_dir  <- file.path(PCA_BASE, tag)
  out_dir  <- file.path(PLOT_BASE, tag)

  eigenvec <- file.path(pca_dir, "pca.eigenvec")
  eigenval <- file.path(pca_dir, "pca.eigenval")

  if (!file.exists(eigenvec)) {
    message("SKIP (no eigenvec): ", eigenvec)
    next
  }

  message("\n", strrep("━", 60))
  message("Dataset : ", tag)

  pca_obj <- read_plink_pca(eigenvec, eigenval)
  message("  Samples : ", nrow(pca_obj$scores))
  message("  PCs     : ", pca_obj$n_pcs)
  message("  PC1 var : ", sprintf("%.2f%%", pca_obj$pct[1]))

  # Scree plot (one per dataset)
  message("  Plotting scree...")
  plot_scree(pca_obj, tag, out_dir)

  # Biplot(s) — one per PC pair
  for (pair in cfg$pc_pairs) {
    pcx <- pair[1]; pcy <- pair[2]

    # Validate PC axes exist in eigenvec
    if (!pcx %in% colnames(pca_obj$scores) || !pcy %in% colnames(pca_obj$scores)) {
      message("  SKIP (PC axes not available): ", pcx, " vs ", pcy)
      next
    }

    message("  Plotting ", pcx, " vs ", pcy, "...")
    plot_pca_biplot(
      pca_obj  = pca_obj,
      cfg      = cfg,
      meta     = meta,
      size_aux = size_aux,
      pcx      = pcx,
      pcy      = pcy,
      out_dir  = out_dir
    )
  }
}

message("\n", strrep("━", 60))
message("All PCA plots complete.")
message("Output directory: ", PLOT_BASE)
