# ============================================================
# INTEGRATED FULL SCRIPT
# Longitudinal scRNA delta-log workflow + trajectories + Reactome ORA
#
# INPUT FILES (all in one folder):
#   C.A_vs_6M.<celltype>.tsv
#   NC.A_vs_6M.<celltype>.tsv
#   C.A_vs_6M.<celltype>.sig_in_NC.tsv
#   NC.A_vs_6M.<celltype>.sig_in_C.tsv
#
# REQUIRED COLUMNS in each file:
#   gene
#   avg_log2FC
#   p_val_adj                 (optional)
#
# ADDITIONAL COLUMNS needed for trajectory plots:
#   avg_expr_group_1
#   avg_expr_group_2
#
# OUTPUT:
#   1) per-cell-type merged delta-log tables
#   2) per-cell-type old-style tables
#   3) Excel workbook with one tab per cell type
#   4) overview delta log all.xlsx
#   5) gene_overlap_summary.xlsx
#   6) scatterplots
#   7) stacked barplot of genes with |delta_log| > 1
#   8) trajectory tables/plots for top 7 cell types + mDCs
#   9) Reactome ORA tables and panels
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(writexl)
  library(fgsea)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

input_dir  <- "F:/chikv-results/reanalysis/delta log"
output_dir <- file.path(input_dir, "delta_log_outputs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "trajectory_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "trajectory_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "reactome_panels_ora_heatmap"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "reactome_panels_ora_heatmap_collapsed"), showWarnings = FALSE, recursive = TRUE)

# ---- Delta-log / summary settings ----
delta_threshold <- 1
trajectory_top_n <- 20
trajectory_force_include <- c(
  "Non_classical_monocytes",
  "Plasmacytoid_dendritic_cells",
  "Myeloid_dendritic_cells",
  "Switched_memory_B_cells"
)
shared_group_cutoff <- 6
shared_group_label <- paste0("Genes_ge", shared_group_cutoff, "_groups")
shared_group_title <- paste0("Genes ≥ ", shared_group_cutoff, " groups")

# ---- ORA settings ----
minSize_ora    <- 5
maxSize_ora    <- 500
ora_fdr_cutoff <- 0.05
bar_n          <- 10

# Remove visually misleading pathogen-specific Reactome terms
exclude_reactome_terms <- c(
  "SARS_COV",
  "COVID",
  "CORONAVIRUS",
  "LEISHMANIA",
  "HIV",
  "INFLUENZA",
  "MEASLES",
  "MALARIA",
  "HEPATITIS",
  "TUBERCULOSIS",
  "PATHOGEN",
  "NERVOUS_SYSTEM"
)

# ---- ORA collapse settings ----
collapse_enabled           <- TRUE
collapse_overlap_threshold <- 0.60
collapse_keep_top_n_pre    <- 40
collapse_bar_n             <- 10

# ---- External ORA resources ----
gmt_file      <- "F:/chikv-results/reactome/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
universe_file <- "F:/chikv-results/reactome/de_genes.txt"

# ---- Plot appearance ORA ----
WRAP_WIDTH_PATH   <- 34
MAX_GENES_PER_CT  <- 10
GENE_LABEL_ANGLE  <- 55

PANEL_W_IN        <- 16.0
PANEL_H_IN        <- 9.0
TITLE_SIZE_PT     <- 19

BAR_BASE_SIZE     <- 14
HEAT_BASE_SIZE    <- 12
PATHWAY_TEXT_SIZE <- 12
GENE_TEXT_SIZE    <- 11
GROUP_TEXT_SIZE   <- 11
LEGEND_TITLE_SIZE <- 13
LEGEND_TEXT_SIZE  <- 12

# ============================================================
# HELPERS
# ============================================================

get_celltype <- function(x) {
  x |>
    basename() |>
    str_remove("^C\\.A_vs_6M\\.") |>
    str_remove("^NC\\.A_vs_6M\\.") |>
    str_remove("\\.sig_in_NC\\.tsv$") |>
    str_remove("\\.sig_in_C\\.tsv$") |>
    str_remove("\\.tsv$")
}

read_deg <- function(path, label) {
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }
  
  df <- read_tsv(path, show_col_types = FALSE)
  
  if (!"gene" %in% colnames(df)) stop("Missing gene column in ", basename(path))
  if (!"avg_log2FC" %in% colnames(df)) stop("Missing avg_log2FC in ", basename(path))
  
  if (!"p_val_adj" %in% colnames(df)) {
    df$p_val_adj <- NA_real_
  }
  
  df %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    rename(
      !!paste0(label, "_log2FC") := avg_log2FC,
      !!paste0(label, "_padj")   := p_val_adj
    )
}

read_deg_expr_full <- function(path, prefix) {
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }
  
  df <- read_tsv(path, show_col_types = FALSE)
  
  needed <- c("gene", "avg_expr_group_1", "avg_expr_group_2", "avg_log2FC")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  }
  
  if (!"p_val_adj" %in% colnames(df)) {
    df$p_val_adj <- NA_real_
  }
  
  df %>%
    select(gene, avg_expr_group_1, avg_expr_group_2, avg_log2FC, p_val_adj) %>%
    rename(
      !!paste0(prefix, "_A_expr")  := avg_expr_group_1,
      !!paste0(prefix, "_M6_expr") := avg_expr_group_2,
      !!paste0(prefix, "_log2FC")  := avg_log2FC,
      !!paste0(prefix, "_padj")    := p_val_adj
    )
}

pad_rows <- function(df, n) {
  if (nrow(df) < n) {
    add_n <- n - nrow(df)
    df[(nrow(df) + 1):n, ] <- NA
  }
  df
}

pretty_celltype <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace("^Total$", "PBMC") %>%
    str_replace("^total cells$", "PBMC") %>%
    str_replace("^Naive CD8 T cells$", "Naive CD8 T") %>%
    str_replace("^Natural killer cells$", "NK cells") %>%
    str_replace("^Intermediate monocytes$", "Intermediate monocytes") %>%
    str_replace("^Classical monocytes$", "Classical monocytes") %>%
    str_replace("^Non classical monocytes$", "Non-classical monocytes") %>%
    str_replace("^Non Vd2 gd T cells$", "Non-Vδ2 γδ T") %>%
    str_replace("^Plasmacytoid dendritic cells$", "pDCs") %>%
    str_replace("^Naive B cells$", "Naive B") %>%
    str_replace("^Switched memory B cells$", "Switched memory B") %>%
    str_replace("^Th2 cells$", "Th2") %>%
    str_replace("^MAIT cells$", "MAIT") %>%
    str_replace("^Naive CD4 T cells$", "Naive CD4 T") %>%
    str_replace("^Myeloid dendritic cells$", "mDCs") %>%
    str_replace("^Central memory CD8 T cells$", "CD8 CM") %>%
    str_replace("^Th1 cells$", "Th1") %>%
    str_replace("^Vd2 gd T cells$", "Vδ2 γδ T") %>%
    str_replace("^T regulatory cells$", "Treg") %>%
    str_replace("^Th1 Th17 cells$", "Th1/Th17") %>%
    str_replace("^Exhausted B cells$", "Exhausted B") %>%
    str_replace("^Th17 cells$", "Th17") %>%
    str_replace("^Non switched memory B cells$", "Non-switched memory B") %>%
    str_replace("^Effector memory CD8 T cells$", "CD8 EM") %>%
    str_replace("^Follicular helper T cells$", "Tfh") %>%
    str_replace("^Terminal effector CD8 T cells$", "CD8 TEMRA") %>%
    str_replace("^Terminal effector CD4 T cells$", "CD4 TEMRA")
}

pretty_pathway <- function(x) {
  x |>
    gsub("^REACTOME[:_ ]*", "", x = _) |>
    gsub("_", " ", x = _) |>
    stringr::str_squish()
}

pretty_ct_ora <- function(ct) {
  out <- ct
  out[out == "total cells"] <- "PBMC"
  out[out == "Terminal effector CD8 T cells"] <- "CD8 TEMRA"
  out[out == "Myeloid_dendritic_cells"] <- "mDCs"
  out[out == shared_group_label] <- shared_group_title
  out <- pretty_celltype(out)
  out
}

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

safe_sig <- function(p) {
  s <- -log10(p)
  s[!is.finite(s)] <- NA_real_
  s
}

filter_unwanted_pathways <- function(df) {
  if (!nrow(df)) return(df)
  
  pattern <- paste(exclude_reactome_terms, collapse = "|")
  
  df %>%
    filter(
      !str_detect(
        toupper(pathway),
        pattern
      )
    )
}

unify_padj <- function(df) {
  if (!nrow(df)) {
    df$padj_u <- numeric(0)
    return(df)
  }
  if ("padj" %in% names(df)) {
    df$padj_u <- suppressWarnings(as.numeric(df$padj))
  } else if ("p.adjust" %in% names(df)) {
    df$padj_u <- suppressWarnings(as.numeric(df$p.adjust))
  } else {
    df$padj_u <- NA_real_
  }
  df
}

normalize_overlap_genes <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

jaccard_index <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (!length(a) && !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

collapse_ora_pathways <- function(df_ct,
                                  pathway_list,
                                  top_n_pre = 40,
                                  overlap_threshold = 0.60) {
  if (!nrow(df_ct)) return(df_ct)
  
  df_ct <- df_ct %>%
    arrange(padj_u) %>%
    slice_head(n = top_n_pre)
  
  keep_idx <- integer(0)
  
  for (i in seq_len(nrow(df_ct))) {
    this_path <- df_ct$pathway[i]
    this_genes <- pathway_list[[this_path]]
    if (is.null(this_genes) || !length(this_genes)) next
    
    redundant <- FALSE
    
    if (length(keep_idx)) {
      for (j in keep_idx) {
        prev_path <- df_ct$pathway[j]
        prev_genes <- pathway_list[[prev_path]]
        if (is.null(prev_genes) || !length(prev_genes)) next
        
        jac <- jaccard_index(this_genes, prev_genes)
        overlap_small <- length(intersect(this_genes, prev_genes)) /
          min(length(unique(this_genes)), length(unique(prev_genes)))
        
        if (jac >= overlap_threshold || overlap_small >= 0.80) {
          redundant <- TRUE
          break
        }
      }
    }
    
    if (!redundant) keep_idx <- c(keep_idx, i)
  }
  
  df_ct[keep_idx, , drop = FALSE]
}

# ============================================================
# DETECT CELL TYPES
# ============================================================

all_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

celltypes <- all_files |>
  basename() |>
  keep(~ str_detect(.x, "^(C|NC)\\.A_vs_6M\\.")) |>
  map_chr(get_celltype) |>
  unique() |>
  sort()

cat("Detected cell types:\n")
print(celltypes)

# ============================================================
# CONTAINERS
# ============================================================

all_merged   <- list()
excel_sheets <- list()
scatter_plots <- list()

# ============================================================
# MAIN LOOP: BUILD DELTA-LOG TABLES
# ============================================================

for (ct in celltypes) {
  
  message("Processing: ", ct)
  
  f_C_main      <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))
  f_NC_main     <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))
  f_C_sig_in_NC <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))
  f_NC_sig_in_C <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv"))
  
  C_main   <- read_deg(f_C_main, "C")
  NC_main  <- read_deg(f_NC_main, "NC")
  C_other  <- read_deg(f_C_sig_in_NC, "C")
  NC_other <- read_deg(f_NC_sig_in_C, "NC")
  
  if (is.null(C_main) || is.null(NC_main) || is.null(C_other) || is.null(NC_other)) {
    next
  }
  
  genes_C  <- C_main$gene
  genes_NC <- NC_main$gene
  
  merged <- full_join(
    C_main %>% rename(C_log2FC_main = C_log2FC, C_padj_main = C_padj),
    C_other %>% rename(C_log2FC_other = C_log2FC, C_padj_other = C_padj),
    by = "gene"
  ) %>%
    full_join(
      NC_main %>% rename(NC_log2FC_main = NC_log2FC, NC_padj_main = NC_padj),
      by = "gene"
    ) %>%
    full_join(
      NC_other %>% rename(NC_log2FC_other = NC_log2FC, NC_padj_other = NC_padj),
      by = "gene"
    ) %>%
    mutate(
      C_log2FC  = coalesce(C_log2FC_main,  C_log2FC_other),
      NC_log2FC = coalesce(NC_log2FC_main, NC_log2FC_other),
      C_padj    = coalesce(C_padj_main,    C_padj_other),
      NC_padj   = coalesce(NC_padj_main,   NC_padj_other),
      sig_in_C  = gene %in% genes_C,
      sig_in_NC = gene %in% genes_NC,
      significance_group = case_when(
        sig_in_C & sig_in_NC  ~ "Both",
        sig_in_C & !sig_in_NC ~ "Chronic only",
        !sig_in_C & sig_in_NC ~ "Non-chronic only",
        TRUE                  ~ "Unexpected"
      ),
      delta_log = C_log2FC - NC_log2FC,
      cell_type = ct
    ) %>%
    select(
      cell_type, gene,
      C_log2FC, NC_log2FC, delta_log,
      C_padj, NC_padj,
      sig_in_C, sig_in_NC,
      significance_group
    ) %>%
    filter(!is.na(C_log2FC) & !is.na(NC_log2FC))
  
  merged_out <- merged %>%
    arrange(desc(delta_log))
  
  write_tsv(
    merged_out,
    file.path(output_dir, "tables", paste0(ct, "_delta_log_merged.tsv"))
  )
  
  # ---- old style table ----
  chronic_tbl <- merged %>%
    filter(sig_in_C) %>%
    arrange(desc(abs(C_log2FC))) %>%
    transmute(
      Chronic = NA_character_,
      gene = gene,
      avg_log2FC = C_log2FC
    )
  
  nc_tbl <- merged %>%
    filter(sig_in_NC) %>%
    arrange(desc(abs(NC_log2FC))) %>%
    transmute(
      `Non Chronic` = NA_character_,
      gene = gene,
      avg_log2FC = NC_log2FC
    )
  
  delta_tbl <- merged %>%
    arrange(desc(delta_log)) %>%
    transmute(
      Column1 = gene,
      `C 6m vs A avg_log2FC`  = C_log2FC,
      `NC 6m vs A avg_log2FC` = NC_log2FC,
      `delta log`             = delta_log
    )
  
  max_n <- max(nrow(chronic_tbl), nrow(nc_tbl), nrow(delta_tbl))
  
  chronic_tbl <- pad_rows(chronic_tbl, max_n)
  nc_tbl      <- pad_rows(nc_tbl, max_n)
  delta_tbl   <- pad_rows(delta_tbl, max_n)
  
  old_style <- bind_cols(
    chronic_tbl,
    nc_tbl,
    tibble(blank = rep(NA_character_, max_n)),
    delta_tbl
  )
  
  colnames(old_style) <- c(
    "Chronic", "gene", "avg_log2FC",
    "Non Chronic", "gene", "avg_log2FC",
    "",
    "Column1",
    "C 6m vs A avg_log2FC",
    "NC 6m vs A avg_log2FC",
    "delta log"
  )
  
  write_tsv(
    old_style,
    file.path(output_dir, "tables", paste0(ct, "_old_style_table.tsv"))
  )
  
  sheet_name <- str_sub(ct, 1, 31)
  excel_sheets[[sheet_name]] <- old_style
  
  # ---- scatterplot ----
  p <- ggplot(
    merged,
    aes(x = NC_log2FC, y = C_log2FC, colour = significance_group)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65") +
    geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 1, colour = "#d73027", linewidth = 0.6, linetype = "dashed") +
    geom_abline(slope = 1, intercept = -1, colour = "#d73027", linewidth = 0.6, linetype = "dashed") +
    geom_point(size = 2, alpha = 0.9) +
    scale_colour_manual(values = c(
      "Both" = "#7b2cbf",
      "Chronic only" = "#f4a261",
      "Non-chronic only" = "#e76f51",
      "Unexpected" = "grey70"
    )) +
    labs(
      title = pretty_celltype(ct),
      x = "NC log2FC (A vs 6M)",
      y = "C log2FC (A vs 6M)",
      colour = NULL
    ) +
    coord_equal(
      xlim = c(-6, 6),
      ylim = c(-6, 6),
      expand = FALSE
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  scatter_plots[[ct]] <- p
  
  ggsave(
    file.path(output_dir, "plots", paste0(ct, "_scatter.png")),
    p,
    width = 7,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  all_merged[[ct]] <- merged_out
}

# ============================================================
# COMBINE ALL CELL TYPES
# ============================================================

combined <- bind_rows(all_merged)

write_tsv(
  combined,
  file.path(output_dir, "all_celltypes_delta_log_merged.tsv")
)

excel_sheets[["all_celltypes_merged"]] <- combined

write_xlsx(
  excel_sheets,
  path = file.path(output_dir, "delta_log_old_style_tables.xlsx")
)

# ============================================================
# OVERVIEW DELTA LOG ALL
# ============================================================

overview_list <- list()
cts <- unique(combined$cell_type)

for (ct in cts) {
  tmp <- combined %>%
    filter(cell_type == ct) %>%
    arrange(desc(delta_log)) %>%
    select(gene, delta_log)
  
  colnames(tmp) <- c(ct, paste0(ct, "_delta_log"))
  overview_list[[ct]] <- tmp
}

max_rows <- max(sapply(overview_list, nrow))

overview_wide <- bind_cols(
  lapply(overview_list, pad_rows, n = max_rows)
)

write_xlsx(
  list("overview delta log all" = overview_wide),
  path = file.path(output_dir, "overview delta log all.xlsx")
)

# ============================================================
# GENE OVERLAP SUMMARY
# ============================================================

filtered <- combined %>%
  filter(!is.na(delta_log), abs(delta_log) > delta_threshold)

counts_per_celltype <- filtered %>%
  group_by(cell_type) %>%
  summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
  arrange(desc(n_genes))

presence_matrix <- filtered %>%
  distinct(gene, cell_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = value,
    values_fill = 0
  ) %>%
  arrange(gene)

values_wide <- filtered %>%
  distinct(gene, cell_type, delta_log) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = delta_log
  ) %>%
  arrange(gene)

presence_counts <- presence_matrix %>%
  mutate(n_celltypes = rowSums(across(-gene)))

values_overlap_only <- values_wide %>%
  inner_join(
    presence_counts %>%
      filter(n_celltypes >= 2) %>%
      select(gene),
    by = "gene"
  )

top_genes_by_presence <- presence_counts %>%
  arrange(desc(n_celltypes), gene)

gene_lists <- lapply(cts, function(ct) {
  filtered %>%
    filter(cell_type == ct) %>%
    arrange(desc(delta_log)) %>%
    pull(gene)
})

names(gene_lists) <- cts
max_len <- max(lengths(gene_lists))

gene_by_celltype_wide <- bind_cols(
  lapply(gene_lists, function(x) {
    length(x) <- max_len
    tibble(value = x)
  })
)

colnames(gene_by_celltype_wide) <- cts

threshold_summary <- tibble(
  parameter = c(
    "threshold",
    "n_rows_combined",
    "n_rows_filtered",
    "n_unique_genes",
    "n_celltypes"
  ),
  value = c(
    delta_threshold,
    nrow(combined),
    nrow(filtered),
    n_distinct(filtered$gene),
    n_distinct(filtered$cell_type)
  )
)

pairs <- combn(cts, 2, simplify = FALSE)

pairwise_overlap <- map_dfr(pairs, function(x) {
  a <- x[1]
  b <- x[2]
  
  genes_a <- filtered %>% filter(cell_type == a) %>% pull(gene) %>% unique()
  genes_b <- filtered %>% filter(cell_type == b) %>% pull(gene) %>% unique()
  
  ov <- length(intersect(genes_a, genes_b))
  un <- length(union(genes_a, genes_b))
  
  tibble(
    CellType_A = a,
    CellType_B = b,
    Overlap = ov,
    Jaccard = ov / un
  )
})

value_correlations <- map_dfr(pairs, function(x) {
  a <- x[1]
  b <- x[2]
  
  tmp <- filtered %>%
    filter(cell_type %in% c(a, b)) %>%
    select(gene, cell_type, delta_log) %>%
    pivot_wider(
      names_from = cell_type,
      values_from = delta_log
    )
  
  tmp2 <- tmp %>% drop_na()
  cor_val <- if (nrow(tmp2) >= 2) cor(tmp2[[2]], tmp2[[3]]) else NA_real_
  
  tibble(
    CellType_A = a,
    CellType_B = b,
    n_overlap = nrow(tmp2),
    correlation = cor_val
  )
})

write_xlsx(
  list(
    counts_per_celltype   = counts_per_celltype,
    presence_matrix       = presence_matrix,
    pairwise_overlap      = pairwise_overlap,
    values_wide           = values_wide,
    values_overlap_only   = values_overlap_only,
    value_correlations    = value_correlations,
    top_genes_by_presence = top_genes_by_presence,
    gene_by_celltype_wide = gene_by_celltype_wide,
    threshold_summary     = threshold_summary
  ),
  path = file.path(output_dir, "gene_overlap_summary.xlsx")
)

# ============================================================
# STACKED BARPLOT OF GENES WITH |delta_log| > 1
# ============================================================

bar_df <- combined %>%
  filter(!is.na(delta_log)) %>%
  mutate(abs_delta_log = abs(delta_log)) %>%
  filter(!cell_type %in% c("total cells", "total_cells", "Total")) %>%
  filter(abs_delta_log > 1) %>%
  mutate(
    delta_bin = case_when(
      abs_delta_log > 1 & abs_delta_log <= 2 ~ "1-2",
      abs_delta_log > 2 & abs_delta_log <= 3 ~ "2-3",
      abs_delta_log > 3 & abs_delta_log <= 4 ~ "3-4",
      abs_delta_log > 4 & abs_delta_log <= 5 ~ "4-5",
      abs_delta_log > 5                      ~ "5+"
    )
  ) %>%
  count(cell_type, delta_bin, name = "n_genes") %>%
  complete(
    cell_type,
    delta_bin = c("1-2", "2-3", "3-4", "4-5", "5+"),
    fill = list(n_genes = 0)
  )

celltype_order <- bar_df %>%
  group_by(cell_type) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(cell_type)

pretty_levels <- pretty_celltype(celltype_order)

bar_df <- bar_df %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    cell_type_pretty = factor(
      pretty_celltype(as.character(cell_type)),
      levels = pretty_levels
    ),
    delta_bin = factor(
      delta_bin,
      levels = c("1-2", "2-3", "3-4", "4-5", "5+")
    )
  )

p_bar_summary <- ggplot(
  bar_df,
  aes(x = cell_type_pretty, y = n_genes, fill = delta_bin)
) +
  geom_col(
    width = 0.82,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = c(
      "1-2" = "#66c2a5",
      "2-3" = "#8da0cb",
      "3-4" = "#fc8d62",
      "4-5" = "#e78ac3",
      "5+"  = "#d73027"
    ),
    breaks = c("5+", "4-5", "3-4", "2-3", "1-2")
  ) +
  labs(
    title = expression("Genes with " * "|" * Delta * "log" * "|" > 1 * " per cell type"),
    x = NULL,
    y = expression("Genes with " * "|" * Delta * "log" * "|" > 1),
    fill = expression("|" * Delta * "log" * "|")
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = file.path(output_dir, "plots", "celltype_delta_log_bin_counts.png"),
  plot = p_bar_summary,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

write_tsv(
  bar_df,
  file.path(output_dir, "tables", "celltype_delta_log_bin_counts.tsv")
)

# ============================================================
# TRAJECTORY PLOTS FOR TOP 7 CELL TYPES + mDCs
# ============================================================

trajectory_top_celltypes <- combined %>%
  filter(!is.na(delta_log), abs(delta_log) > 1) %>%
  count(cell_type, name = "n_genes_abs_delta_gt1") %>%
  arrange(desc(n_genes_abs_delta_gt1)) %>%
  slice_head(n = 7)

write_tsv(
  trajectory_top_celltypes,
  file.path(output_dir, "tables", "top7_celltypes_by_abs_delta_gt1.tsv")
)

trajectory_celltypes_auto <- unique(c(
  trajectory_top_celltypes$cell_type,
  trajectory_force_include
))

trajectory_excel <- list()

for (ct in trajectory_celltypes_auto) {
  
  message("Creating trajectory plots for: ", ct)
  
  merged_out_ct <- combined %>%
    filter(cell_type == ct) %>%
    arrange(desc(delta_log))
  
  if (nrow(merged_out_ct) == 0) {
    warning("No merged data found for: ", ct)
    next
  }
  
  top_genes_tbl <- merged_out_ct %>%
    mutate(abs_delta_log = abs(delta_log)) %>%
    arrange(desc(abs_delta_log), desc(delta_log)) %>%
    slice_head(n = trajectory_top_n) %>%
    select(gene, delta_log, abs_delta_log)
  
  f_C_main      <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))
  f_NC_main     <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))
  f_C_sig_in_NC <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))
  f_NC_sig_in_C <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv"))
  
  C_main_expr   <- read_deg_expr_full(f_C_main, "C")
  NC_main_expr  <- read_deg_expr_full(f_NC_main, "NC")
  C_other_expr  <- read_deg_expr_full(f_C_sig_in_NC, "C")
  NC_other_expr <- read_deg_expr_full(f_NC_sig_in_C, "NC")
  
  if (is.null(C_main_expr) || is.null(NC_main_expr) || is.null(C_other_expr) || is.null(NC_other_expr)) {
    next
  }
  
  expr_merged <- full_join(
    C_main_expr %>% rename(
      C_A_expr_main  = C_A_expr,
      C_M6_expr_main = C_M6_expr
    ),
    C_other_expr %>% rename(
      C_A_expr_other  = C_A_expr,
      C_M6_expr_other = C_M6_expr
    ),
    by = "gene"
  ) %>%
    full_join(
      NC_main_expr %>% rename(
        NC_A_expr_main  = NC_A_expr,
        NC_M6_expr_main = NC_M6_expr
      ),
      by = "gene"
    ) %>%
    full_join(
      NC_other_expr %>% rename(
        NC_A_expr_other  = NC_A_expr,
        NC_M6_expr_other = NC_M6_expr
      ),
      by = "gene"
    ) %>%
    mutate(
      C_A_expr   = coalesce(C_A_expr_main,  C_A_expr_other),
      C_M6_expr  = coalesce(C_M6_expr_main, C_M6_expr_other),
      NC_A_expr  = coalesce(NC_A_expr_main,  NC_A_expr_other),
      NC_M6_expr = coalesce(NC_M6_expr_main, NC_M6_expr_other)
    ) %>%
    select(gene, C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr)
  
  top_genes_with_expr <- top_genes_tbl %>%
    left_join(expr_merged, by = "gene") %>%
    mutate(
      has_complete_expr = !is.na(C_A_expr) & !is.na(C_M6_expr) &
        !is.na(NC_A_expr) & !is.na(NC_M6_expr)
    )
  
  write_tsv(
    top_genes_with_expr,
    file.path(output_dir, "trajectory_tables", paste0(ct, "_top", trajectory_top_n, "_genes_from_main_delta_log.tsv"))
  )
  
  trajectory_excel[[str_sub(ct, 1, 31)]] <- top_genes_with_expr
  
  plot_tbl <- top_genes_with_expr %>%
    filter(has_complete_expr)
  
  if (nrow(plot_tbl) == 0) next
  
  long_df <- plot_tbl %>%
    select(gene, delta_log, abs_delta_log, C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr) %>%
    pivot_longer(
      cols = c(C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr),
      names_to = "group_time",
      values_to = "avg_expr"
    ) %>%
    mutate(
      Outcome = case_when(
        str_starts(group_time, "C_")  ~ "Chronic",
        str_starts(group_time, "NC_") ~ "Non-chronic"
      ),
      Time = case_when(
        str_detect(group_time, "_A_expr$")  ~ "Acute",
        str_detect(group_time, "_M6_expr$") ~ "6M"
      ),
      Time = factor(Time, levels = c("Acute", "6M")),
      gene = factor(gene, levels = plot_tbl$gene)
    )
  
  p_traj <- ggplot(
    long_df,
    aes(x = Time, y = avg_expr, group = Outcome, colour = Outcome)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = c(
      "Chronic" = "#d73027",
      "Non-chronic" = "#4575b4"
    )) +
    labs(
      title = paste0(pretty_celltype(ct), ": top ", trajectory_top_n, " genes by |Δlog|"),
      x = NULL,
      y = "Average expression",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(output_dir, "trajectory_plots", paste0(ct, "_top", trajectory_top_n, "_trajectory_plot.png")),
    plot = p_traj,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
  )
}

if (length(trajectory_excel) > 0) {
  write_xlsx(
    trajectory_excel,
    path = file.path(output_dir, "top7plus_mDC_top20_trajectory_genes.xlsx")
  )
}

# ============================================================
# CUSTOM TRAJECTORY FIGURES FOR SELECTED MAIN-FIGURE GENES
# ============================================================

dir.create(
  file.path(output_dir, "trajectory_custom_selected"),
  recursive = TRUE,
  showWarnings = FALSE
)

make_custom_celltype_trajectory <- function(ct, genes_keep, fig_name, fig_title, ncol = 2) {
  
  message("Creating custom trajectory figure for: ", ct)
  
  f_C_main      <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))
  f_NC_main     <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))
  f_C_sig_in_NC <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))
  f_NC_sig_in_C <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv"))
  
  C_main_expr   <- read_deg_expr_full(f_C_main, "C")
  NC_main_expr  <- read_deg_expr_full(f_NC_main, "NC")
  C_other_expr  <- read_deg_expr_full(f_C_sig_in_NC, "C")
  NC_other_expr <- read_deg_expr_full(f_NC_sig_in_C, "NC")
  
  expr_merged <- full_join(
    C_main_expr %>% rename(C_A_expr_main = C_A_expr, C_M6_expr_main = C_M6_expr, C_log2FC_main = C_log2FC),
    C_other_expr %>% rename(C_A_expr_other = C_A_expr, C_M6_expr_other = C_M6_expr, C_log2FC_other = C_log2FC),
    by = "gene"
  ) %>%
    full_join(
      NC_main_expr %>% rename(NC_A_expr_main = NC_A_expr, NC_M6_expr_main = NC_M6_expr, NC_log2FC_main = NC_log2FC),
      by = "gene"
    ) %>%
    full_join(
      NC_other_expr %>% rename(NC_A_expr_other = NC_A_expr, NC_M6_expr_other = NC_M6_expr, NC_log2FC_other = NC_log2FC),
      by = "gene"
    ) %>%
    mutate(
      C_A_expr   = coalesce(C_A_expr_main, C_A_expr_other),
      C_M6_expr  = coalesce(C_M6_expr_main, C_M6_expr_other),
      NC_A_expr  = coalesce(NC_A_expr_main, NC_A_expr_other),
      NC_M6_expr = coalesce(NC_M6_expr_main, NC_M6_expr_other),
      C_log2FC   = coalesce(C_log2FC_main, C_log2FC_other),
      NC_log2FC  = coalesce(NC_log2FC_main, NC_log2FC_other),
      delta_log  = C_log2FC - NC_log2FC,
      gene_label = paste0(gene, "\nΔlog=", round(delta_log, 2))
    ) %>%
    select(gene, gene_label, C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr, C_log2FC, NC_log2FC, delta_log) %>%
    filter(gene %in% genes_keep)
  
  label_levels <- expr_merged %>%
    mutate(gene_order = match(gene, genes_keep)) %>%
    arrange(gene_order) %>%
    pull(gene_label)
  
  expr_merged <- expr_merged %>%
    mutate(gene_label = factor(gene_label, levels = label_levels))
  
  write_tsv(
    expr_merged,
    file.path(output_dir, "trajectory_custom_selected", paste0(fig_name, "_data.tsv"))
  )
  
  plot_df <- expr_merged %>%
    pivot_longer(
      cols = c(C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr),
      names_to = "group_time",
      values_to = "avg_expr"
    ) %>%
    mutate(
      Outcome = case_when(
        str_starts(group_time, "C_")  ~ "Chronic",
        str_starts(group_time, "NC_") ~ "Non-chronic"
      ),
      Time = case_when(
        str_detect(group_time, "_A_expr$")  ~ "Acute",
        str_detect(group_time, "_M6_expr$") ~ "6M"
      ),
      Time = factor(Time, levels = c("Acute", "6M")),
      Outcome = factor(Outcome, levels = c("Non-chronic", "Chronic"))
    ) %>%
    filter(!is.na(avg_expr))
  
  p <- ggplot(
    plot_df,
    aes(x = Time, y = avg_expr, group = Outcome, colour = Outcome)
  ) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 2.4) +
    facet_wrap(~ gene_label, scales = "free_y", ncol = ncol) +
    scale_colour_manual(values = c(
      "Non-chronic" = "#4575b4",
      "Chronic" = "#d73027"
    )) +
    labs(
      title = fig_title,
      x = NULL,
      y = "Average expression",
      colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(output_dir, "trajectory_custom_selected", paste0(fig_name, ".png")),
    plot = p,
    width = 7,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(output_dir, "trajectory_custom_selected", paste0(fig_name, ".pdf")),
    plot = p,
    width = 7,
    height = 7,
    bg = "white"
  )
  
  return(p)
}
# ---- Non-classical monocytes selected genes ----
nonclassical_selected_genes <- c(
  "SMAD3",
  "HLA-DPA1",
  "FCGR3A",
  "RHOA",
  "SH3RF3",
  "HIPK2"
)

p_traj_nonclassical <- make_custom_celltype_trajectory(
  ct = "Non_classical_monocytes",
  genes_keep = nonclassical_selected_genes,
  fig_name = "Non_classical_monocytes_selected_gene_trajectories",
  fig_title = "Non-classical monocytes: selected divergent trajectories",
  ncol = 2
)

# ---- pDC selected genes ----
pdc_selected_genes <- c(
  "IFI44L",
  "TCF4",
  "LILRA4",
  "GZMB",
  "NFKB2",
  "HLA-DPA1"
)

p_traj_pdc <- make_custom_celltype_trajectory(
  ct = "Plasmacytoid_dendritic_cells",
  genes_keep = pdc_selected_genes,
  fig_name = "pDC_selected_gene_trajectories",
  fig_title = "pDCs: selected divergent trajectories",
  ncol = 2
)
# ============================================================
# REACTOME ORA
# selected trajectory cell types + genes present in >= 4 non-PBMC groups
# ============================================================

message("\n=== Reactome ORA setup ===")

gmt_lines <- readLines(gmt_file)
gmt_split <- strsplit(gmt_lines, "\t")
gs_react <- lapply(gmt_split, function(x) unique(toupper(x[-c(1, 2)])))
names(gs_react) <- vapply(gmt_split, `[`, character(1), 1)

tested_gene_universe <- readLines(universe_file)
tested_gene_universe <- tested_gene_universe |>
  stringr::str_squish() |>
  toupper()
tested_gene_universe <- unique(tested_gene_universe[tested_gene_universe != ""])

gs_react <- lapply(gs_react, function(g) intersect(g, tested_gene_universe))
gs_react <- gs_react[lengths(gs_react) >= minSize_ora]

message("Pathways loaded after universe filtering: ", length(gs_react))
message("Universe size: ", length(tested_gene_universe))

long_ora <- combined %>%
  transmute(
    CellType = cell_type,
    Gene = toupper(gene),
    delta_log = delta_log
  ) %>%
  distinct(CellType, Gene, .keep_all = TRUE)

presence_tbl <- long_ora %>%
  filter(
    !CellType %in% c("total cells", "total_cells","Total"),
    !is.na(delta_log),
    abs(delta_log) >= delta_threshold
  ) %>%
  distinct(CellType, Gene)

genes_ge4 <- presence_tbl %>%
  count(Gene, name = "n_groups") %>%
  filter(n_groups >= shared_group_cutoff) %>%
  pull(Gene) %>%
  unique()

message(
  "Genes with |delta_log| >= ", delta_threshold,
  " in >= ", shared_group_cutoff,
  " non-PBMC groups: ", length(genes_ge4)
)

# ============================================================
# AVERAGE EXPRESSION FOR GENES PRESENT IN >=4 CELL GROUPS
# ============================================================

message("Building average-expression table for genes present in >=4 cell groups...")

shared_ge4_expr_list <- list()

for (ct in unique(combined$cell_type)) {
  
  f_C_main      <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))
  f_NC_main     <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))
  f_C_sig_in_NC <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))
  f_NC_sig_in_C <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv"))
  
  C_main_expr   <- read_deg_expr_full(f_C_main, "C")
  NC_main_expr  <- read_deg_expr_full(f_NC_main, "NC")
  C_other_expr  <- read_deg_expr_full(f_C_sig_in_NC, "C")
  NC_other_expr <- read_deg_expr_full(f_NC_sig_in_C, "NC")
  
  if (is.null(C_main_expr) || is.null(NC_main_expr) ||
      is.null(C_other_expr) || is.null(NC_other_expr)) {
    next
  }
  
  expr_merged_ct <- full_join(
    C_main_expr %>% rename(
      C_A_expr_main  = C_A_expr,
      C_M6_expr_main = C_M6_expr
    ),
    C_other_expr %>% rename(
      C_A_expr_other  = C_A_expr,
      C_M6_expr_other = C_M6_expr
    ),
    by = "gene"
  ) %>%
    full_join(
      NC_main_expr %>% rename(
        NC_A_expr_main  = NC_A_expr,
        NC_M6_expr_main = NC_M6_expr
      ),
      by = "gene"
    ) %>%
    full_join(
      NC_other_expr %>% rename(
        NC_A_expr_other  = NC_A_expr,
        NC_M6_expr_other = NC_M6_expr
      ),
      by = "gene"
    ) %>%
    mutate(
      C_A_expr   = coalesce(C_A_expr_main,  C_A_expr_other),
      C_M6_expr  = coalesce(C_M6_expr_main, C_M6_expr_other),
      NC_A_expr  = coalesce(NC_A_expr_main,  NC_A_expr_other),
      NC_M6_expr = coalesce(NC_M6_expr_main, NC_M6_expr_other),
      cell_type = ct,
      Gene = toupper(gene)
    ) %>%
    select(
      cell_type,
      Gene,
      C_A_expr,
      C_M6_expr,
      NC_A_expr,
      NC_M6_expr
    )
  
  shared_ge4_expr_list[[ct]] <- expr_merged_ct
}

shared_ge4_expr_long <- bind_rows(shared_ge4_expr_list) %>%
  filter(Gene %in% genes_ge4) %>%
  filter(!cell_type %in% c("total cells", "total_cells", "Total"))%>%
  left_join(
    combined %>%
      transmute(
        cell_type,
        Gene = toupper(gene),
        C_log2FC,
        NC_log2FC,
        delta_log,
        significance_group
      ),
    by = c("cell_type", "Gene")
  ) %>%
  filter(!is.na(delta_log), abs(delta_log) >= delta_threshold) %>%
  arrange(Gene, cell_type)

shared_ge4_expr_wide <- shared_ge4_expr_long %>%
  mutate(cell_type_pretty = pretty_celltype(cell_type)) %>%
  select(
    Gene,
    cell_type_pretty,
    C_A_expr,
    C_M6_expr,
    NC_A_expr,
    NC_M6_expr,
    C_log2FC,
    NC_log2FC,
    delta_log,
    significance_group
  ) %>%
  pivot_wider(
    names_from = cell_type_pretty,
    values_from = c(
      C_A_expr,
      C_M6_expr,
      NC_A_expr,
      NC_M6_expr,
      C_log2FC,
      NC_log2FC,
      delta_log,
      significance_group
    ),
    names_glue = "{cell_type_pretty}_{.value}"
  )

write_tsv(
  shared_ge4_expr_long,
  file.path(output_dir, "tables", "genes_ge4_average_expression_long.tsv")
)

write_xlsx(
  list(
    genes_ge4_expression_long = shared_ge4_expr_long,
    genes_ge4_expression_wide = shared_ge4_expr_wide
  ),
  path = file.path(output_dir, "genes_ge4_average_expression.xlsx")
)
# ============================================================
# CUSTOM TRAJECTORY FIGURE FOR SELECTED GENES PRESENT IN >=6 CELL GROUPS
# ============================================================
shared_selected_genes <- c(
  "RICTOR",
  "MECP2",
  "IFI44L",
  "BIRC3",
  "S100A4",
  "ACTG1",
  "HLA-DPA1",
  "HOPX"
)

shared_selected_long <- shared_ge4_expr_long %>%
  filter(Gene %in% shared_selected_genes) %>%
  mutate(
    cell_type_pretty = pretty_celltype(cell_type),
    
    cell_type_pretty = factor(
      cell_type_pretty,
      levels = c(
        "Classical monocytes",
        "Intermediate monocytes",
        "Non-classical monocytes",
        "pDCs",
        "mDCs",
        "NK cells",
        "MAIT",
        "Vδ2 γδ T",
        "Non-Vδ2 γδ T",
        "Naive CD8 T",
        "CD8 CM",
        "CD8 EM",
        "CD8 TEMRA",
        "Naive CD4 T",
        "Treg",
        "Th1",
        "Th2",
        "Th17",
        "Th1/Th17",
        "Tfh",
        "CD4 TEMRA",
        "Naive B",
        "Switched memory B",
        "Non-switched memory B",
        "Exhausted B"
      )
    ),
    
    Gene = factor(Gene, levels = shared_selected_genes),
    
    Celltype_label = paste0(
      cell_type_pretty,
      "\nΔlog=", round(delta_log, 2)
    )
  ) %>%
  pivot_longer(
    cols = c(C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr),
    names_to = "group_time",
    values_to = "avg_expr"
  ) %>%
  mutate(
    Outcome = case_when(
      str_starts(group_time, "C_")  ~ "Chronic",
      str_starts(group_time, "NC_") ~ "Non-chronic"
    ),
    Time = case_when(
      str_detect(group_time, "_A_expr$")  ~ "Acute",
      str_detect(group_time, "_M6_expr$") ~ "6M"
    ),
    Time = factor(Time, levels = c("Acute", "6M")),
    Outcome = factor(Outcome, levels = c("Non-chronic", "Chronic"))
  ) %>%
  filter(!is.na(avg_expr))

write_tsv(
  shared_selected_long,
  file.path(output_dir, "trajectory_custom_selected", "shared_ge6_selected_gene_trajectories_long.tsv")
)

delta_annot <- shared_selected_long %>%
  distinct(Gene, cell_type_pretty, delta_log) %>%
  mutate(
    label = paste0("Δlog=", round(delta_log, 2)),
    x = 1.5,
    y = Inf
  )

celltype_labeller <- labeller(
  cell_type_pretty = label_wrap_gen(width = 14)
)

p_shared_selected <- ggplot(
  shared_selected_long,
  aes(x = Time, y = avg_expr, group = Outcome, colour = Outcome)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  geom_text(
    data = delta_annot,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 2.6,
    colour = "grey25"
  ) +
  facet_grid(
    Gene ~ cell_type_pretty,
    scales = "free_y",
    labeller = celltype_labeller
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  scale_colour_manual(values = c(
    "Non-chronic" = "#4575b4",
    "Chronic" = "#d73027"
  )) +
  labs(
    title = "Selected conserved divergent trajectories across immune cell types",
    x = NULL,
    y = "Average expression",
    colour = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    strip.text.x = element_text(face = "bold", size = 8),
    strip.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 20, b = 20, l = 20)
  )

ggsave(
  filename = file.path(output_dir, "trajectory_custom_selected", "shared_ge6_selected_gene_trajectories.png"),
  plot = p_shared_selected,
  width = 18,
  height = 10,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "trajectory_custom_selected", "shared_ge6_selected_gene_trajectories.pdf"),
  plot = p_shared_selected,
  width = 18,
  height = 10,
  bg = "white",
  limitsize = FALSE
)
# ============================================================
# HEATMAP: GENES PRESENT IN >=4 CELL GROUPS
# delta_log across immune cell types
# ============================================================

dir.create(
  file.path(output_dir, "shared_gene_heatmaps"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Matrix: rows = genes, columns = cell types, values = delta_log
shared_heatmap_mat <- shared_ge4_expr_long %>%
  select(Gene, cell_type, delta_log) %>%
  distinct(Gene, cell_type, .keep_all = TRUE) %>%
  mutate(cell_type_pretty = pretty_celltype(cell_type)) %>%
  select(Gene, cell_type_pretty, delta_log) %>%
  pivot_wider(
    names_from = cell_type_pretty,
    values_from = delta_log,
    values_fill = 0
  )


# Convert to matrix-like data frame
shared_heatmap_df <- as.data.frame(shared_heatmap_mat)
rownames(shared_heatmap_df) <- shared_heatmap_df$Gene
shared_heatmap_df$Gene <- NULL

# remove PBMC / total cells if present
if ("total cells" %in% colnames(shared_heatmap_df)) {
  shared_heatmap_df <- shared_heatmap_df[, colnames(shared_heatmap_df) != "total cells", drop = FALSE]
}
if ("PBMC" %in% colnames(shared_heatmap_df)) {
  shared_heatmap_df <- shared_heatmap_df[, colnames(shared_heatmap_df) != "PBMC", drop = FALSE]
}

if ("Total" %in% colnames(shared_heatmap_df)) {
  shared_heatmap_df <- shared_heatmap_df[, colnames(shared_heatmap_df) != "Total", drop = FALSE]
}
# Keep genes present in >=4 cell types after filtering
shared_heatmap_df <- shared_heatmap_df[
  rowSums(shared_heatmap_df != 0, na.rm = TRUE) >= 4,
  ,
  drop = FALSE
]

# Order genes by presence, then mean absolute delta_log
gene_presence <- rowSums(shared_heatmap_df != 0, na.rm = TRUE)
gene_mean_abs <- rowMeans(abs(shared_heatmap_df), na.rm = TRUE)

gene_order <- order(
  gene_presence,
  gene_mean_abs,
  decreasing = TRUE
)

shared_heatmap_df <- shared_heatmap_df[gene_order, , drop = FALSE]

# Biological column order, using pretty names
desired_order <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "Non-classical monocytes",
  "pDCs",
  "mDCs",
  "NK cells",
  "MAIT",
  "Vδ2 γδ T",      
  "Non-Vδ2 γδ T",
  "Naive CD8 T",
  "CD8 CM",
  "CD8 EM",
  "CD8 TEMRA",
  "Naive CD4 T",
  "Treg",
  "Th1",
  "Th2",
  "Th17",
  "Th1/Th17",
  "Tfh",
  "CD4 TEMRA",
  "Naive B",
  "Switched memory B",
  "Non-switched memory B",
  "Exhausted B"
)

desired_order <- desired_order[desired_order %in% colnames(shared_heatmap_df)]

# keep any columns not listed at the end
remaining_cols <- setdiff(colnames(shared_heatmap_df), desired_order)
final_col_order <- c(desired_order, remaining_cols)

shared_heatmap_df <- shared_heatmap_df[, final_col_order, drop = FALSE]

# Long format
shared_heatmap_df$Gene <- rownames(shared_heatmap_df)

shared_heatmap_long <- shared_heatmap_df %>%
  as_tibble() %>%
  pivot_longer(
    cols = -Gene,
    names_to = "CellType",
    values_to = "delta_log"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = rev(rownames(shared_heatmap_df))),
    CellType = factor(CellType, levels = final_col_order)
  )

max_abs <- max(abs(shared_heatmap_long$delta_log), na.rm = TRUE)

p_shared_heatmap <- ggplot(
  shared_heatmap_long,
  aes(x = CellType, y = Gene, fill = delta_log)
) +
  geom_tile(color = "grey92", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "delta log2FC"
  ) +
  labs(
    title = "Genes with conserved trajectory divergence across immune cell types"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 23),
    axis.text.y = element_text(size = 16),
    axis.title = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 20, b = 80, l = 80),
    legend.position = "right",
    panel.grid = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.key.height = unit(1.8, "cm")
  )

ggsave(
  filename = file.path(output_dir, "shared_gene_heatmaps", "genes_ge4_delta_log_heatmap.png"),
  plot = p_shared_heatmap,
  width = 18,
  height = max(8, 0.32 * nrow(shared_heatmap_df)),
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "shared_gene_heatmaps", "genes_ge4_delta_log_heatmap.pdf"),
  plot = p_shared_heatmap,
  width = 18,
  height = max(8, 0.32 * nrow(shared_heatmap_df)),
  bg = "white",
  limitsize = FALSE
)

write_tsv(
  shared_heatmap_long,
  file.path(output_dir, "shared_gene_heatmaps", "genes_ge4_delta_log_heatmap_long.tsv")
)

write_xlsx(
  list(
    heatmap_matrix = shared_heatmap_df,
    heatmap_long = shared_heatmap_long
  ),
  path = file.path(output_dir, "shared_gene_heatmaps", "genes_ge4_delta_log_heatmap_data.xlsx")
)

# ============================================================
# TRAJECTORY FIGURES FOR TOP 25 GENES PRESENT IN >=4 CELL GROUPS
# one figure per gene, showing all cell types where present
# ============================================================

dir.create(
  file.path(output_dir, "trajectory_shared_ge4_genes"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Rank shared genes:
# priority = present in most cell types, then strongest max |delta_log|,
# then strongest mean |delta_log|
top25_ge4_genes <- shared_ge4_expr_long %>%
  group_by(Gene) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),
    max_abs_delta_log = max(abs(delta_log), na.rm = TRUE),
    mean_abs_delta_log = mean(abs(delta_log), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_celltypes),
    desc(max_abs_delta_log),
    desc(mean_abs_delta_log),
    Gene
  ) %>%
  slice_head(n = 25)

write_tsv(
  top25_ge4_genes,
  file.path(output_dir, "tables", "top25_genes_ge4_ranked.tsv")
)

top25_shared_expr <- shared_ge4_expr_long %>%
  filter(Gene %in% top25_ge4_genes$Gene) %>%
  left_join(top25_ge4_genes, by = "Gene") %>%
  mutate(
    cell_type_pretty = pretty_celltype(cell_type),
    cell_type_pretty = factor(
      cell_type_pretty,
      levels = pretty_celltype(
        shared_ge4_expr_long %>%
          count(cell_type, name = "n") %>%
          arrange(desc(n)) %>%
          pull(cell_type)
      )
    )
  )

write_xlsx(
  list(
    top25_genes_ranked = top25_ge4_genes,
    top25_expression_long = top25_shared_expr
  ),
  path = file.path(output_dir, "top25_ge4_gene_trajectories_data.xlsx")
)

# long plotting table
top25_shared_long <- top25_shared_expr %>%
  select(
    Gene,
    cell_type,
    cell_type_pretty,
    C_A_expr,
    C_M6_expr,
    NC_A_expr,
    NC_M6_expr,
    C_log2FC,
    NC_log2FC,
    delta_log,
    n_celltypes,
    max_abs_delta_log,
    mean_abs_delta_log
  ) %>%
  pivot_longer(
    cols = c(C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr),
    names_to = "group_time",
    values_to = "avg_expr"
  ) %>%
  mutate(
    Outcome = case_when(
      str_starts(group_time, "C_")  ~ "Chronic",
      str_starts(group_time, "NC_") ~ "Non-chronic"
    ),
    Time = case_when(
      str_detect(group_time, "_A_expr$")  ~ "Acute",
      str_detect(group_time, "_M6_expr$") ~ "6M"
    ),
    Time = factor(Time, levels = c("Acute", "6M")),
    Outcome = factor(Outcome, levels = c("Non-chronic", "Chronic"))
  )

# one figure per gene
for (g in top25_ge4_genes$Gene) {
  
  plot_df <- top25_shared_long %>%
    filter(Gene == g, !is.na(avg_expr))
  
  if (nrow(plot_df) == 0) next
  
  gene_info <- top25_ge4_genes %>%
    filter(Gene == g)
  
  p_gene <- ggplot(
    plot_df,
    aes(x = Time, y = avg_expr, group = Outcome, colour = Outcome)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.1) +
    facet_wrap(~ cell_type_pretty, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = c(
      "Non-chronic" = "#4575b4",
      "Chronic" = "#d73027"
    )) +
    labs(
      title = paste0(
        g,
        " trajectories across cell types"
      ),
      subtitle = paste0(
        "Present in ", gene_info$n_celltypes,
        " cell groups | max |delta log| = ",
        round(gene_info$max_abs_delta_log, 2),
        " | mean |delta log| = ",
        round(gene_info$mean_abs_delta_log, 2)
      ),
      x = NULL,
      y = "Average expression",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(
      output_dir,
      "trajectory_shared_ge4_genes",
      paste0("trajectory_", safe(g), "_ge4_celltypes.png")
    ),
    plot = p_gene,
    width = 12,
    height = max(5, 2.2 * ceiling(length(unique(plot_df$cell_type_pretty)) / 4)),
    dpi = 300,
    bg = "white"
  )
}

run_ora <- function(ct) {
  df_ct <- long_ora |>
    filter(CellType == ct, !is.na(delta_log)) |>
    distinct(Gene, .keep_all = TRUE)
  
  selected <- df_ct |>
    filter(abs(delta_log) >= delta_threshold) |>
    pull(Gene) |>
    unique()
  
  universe <- tested_gene_universe
  selected <- intersect(selected, universe)
  
  message("  ", ct, ": selected=", length(selected), " | universe=", length(universe))
  
  if (length(selected) < 2 || length(universe) < 10) {
    return(tibble())
  }
  
  res <- fgsea::fora(
    pathways = gs_react,
    genes    = selected,
    universe = universe,
    minSize  = minSize_ora,
    maxSize  = maxSize_ora
  ) |>
    as_tibble() |>
    unify_padj()
  
  res <- filter_unwanted_pathways(res)
  
  if (!nrow(res)) return(tibble())
  if (!"overlapGenes" %in% names(res)) res$overlapGenes <- vector("list", nrow(res))
  
  res <- res |>
    mutate(
      overlapGenes = purrr::map(overlapGenes, normalize_overlap_genes)
    ) |>
    filter(!is.na(padj_u)) |>
    filter(padj_u <= ora_fdr_cutoff) |>
    arrange(padj_u) |>
    mutate(
      CellType    = ct,
      selected_n  = length(selected),
      universe_n  = length(universe),
      effect      = log2((overlap / size) / (length(selected) / length(universe))),
      genes_shown = purrr::map_int(overlapGenes, length)
    )
  
  res
}

run_ora_shared_ge4 <- function(shared_genes, universe, pathways) {
  selected <- intersect(unique(shared_genes), universe)
  
  message("  ", shared_group_title, ": selected=", length(selected), " | universe=", length(universe))
  
  if (length(selected) < 2 || length(universe) < 10) {
    return(tibble())
  }
  
  res <- fgsea::fora(
    pathways = pathways,
    genes    = selected,
    universe = universe,
    minSize  = minSize_ora,
    maxSize  = maxSize_ora
  ) %>%
    as_tibble() %>%
    unify_padj()
  
  res <- filter_unwanted_pathways(res)
  
  if (!nrow(res)) return(tibble())
  if (!"overlapGenes" %in% names(res)) res$overlapGenes <- vector("list", nrow(res))
  
  res <- res %>%
    mutate(
      overlapGenes = purrr::map(overlapGenes, normalize_overlap_genes)
    ) %>%
    filter(!is.na(padj_u)) %>%
    filter(padj_u <= ora_fdr_cutoff) %>%
    arrange(padj_u) %>%
    mutate(
      CellType    = shared_group_label,
      selected_n  = length(selected),
      universe_n  = length(universe),
      effect      = log2((overlap / size) / (length(selected) / length(universe))),
      genes_shown = purrr::map_int(overlapGenes, length)
    )
  
  res
}
ora_targets <- unique(c(
  trajectory_celltypes_auto,
  "Non_classical_monocytes",
  "Plasmacytoid_dendritic_cells",
  "Myeloid_dendritic_cells",
  "Switched_memory_B_cells"
))

ora_targets <- intersect(ora_targets, unique(long_ora$CellType))

message("ORA targets are: ", paste(ora_targets, collapse = " | "))
message("Shared ORA label is: ", shared_group_label)

message("\n=== Running ORA ===")
ORA_main <- purrr::map_dfr(ora_targets, run_ora)
ORA_shared_ge4 <- run_ora_shared_ge4(
  shared_genes = genes_ge4,
  universe     = tested_gene_universe,
  pathways     = gs_react
)

ORA_results <- bind_rows(ORA_main, ORA_shared_ge4)

if (!nrow(ORA_results)) {
  warning("No ORA results passed the filters.")
}

if (collapse_enabled && nrow(ORA_results)) {
  ORA_results_collapsed <- ORA_results %>%
    group_split(CellType) %>%
    purrr::map_dfr(function(df_ct) {
      collapse_ora_pathways(
        df_ct = df_ct,
        pathway_list = gs_react,
        top_n_pre = collapse_keep_top_n_pre,
        overlap_threshold = collapse_overlap_threshold
      )
    })
} else {
  ORA_results_collapsed <- ORA_results
}

genes_long_export <- ORA_results |>
  filter(lengths(overlapGenes) > 0) |>
  select(CellType, pathway, effect, padj_u, overlap, size,
         selected_n, universe_n, overlapGenes) |>
  tidyr::unnest_longer(overlapGenes, values_to = "Gene") |>
  mutate(
    Pathway = pretty_pathway(pathway),
    sig     = -log10(padj_u)
  )

results_export <- ORA_results |>
  mutate(overlapGenes = purrr::map_chr(overlapGenes, paste, collapse = "; "))

genes_long_export_collapsed <- ORA_results_collapsed |>
  filter(lengths(overlapGenes) > 0) |>
  select(CellType, pathway, effect, padj_u, overlap, size,
         selected_n, universe_n, overlapGenes) |>
  tidyr::unnest_longer(overlapGenes, values_to = "Gene") |>
  mutate(
    Pathway = pretty_pathway(pathway),
    sig     = -log10(padj_u)
  )

results_export_collapsed <- ORA_results_collapsed |>
  mutate(overlapGenes = purrr::map_chr(overlapGenes, paste, collapse = "; "))

write_xlsx(
  list(
    ORA_results_standard        = results_export,
    ORA_results_collapsed       = results_export_collapsed,
    genes_in_pathways_standard  = genes_long_export,
    genes_in_pathways_collapsed = genes_long_export_collapsed,
    params = tibble(
      delta_cutoff = delta_threshold,
      minSize_ora = minSize_ora,
      maxSize_ora = maxSize_ora,
      ora_fdr_cutoff = ora_fdr_cutoff,
      collapse_enabled = collapse_enabled,
      collapse_overlap_threshold = collapse_overlap_threshold,
      collapse_keep_top_n_pre = collapse_keep_top_n_pre,
      universe_type = "All genes tested in the scRNA analysis (de_genes.txt)",
      targets = paste(c(ora_targets, shared_group_label), collapse = "; ")
    )
  ),
  path = file.path(output_dir, "reactome_ora_results_with_collapsed.xlsx")
)

# ============================================================
# ORA PLOT HELPERS
# ============================================================

logfc_long <- combined %>%
  transmute(
    CellType = cell_type,
    Gene = toupper(gene),
    log2FC_chronic = C_log2FC,
    log2FC_nonchronic = NC_log2FC,
    delta_log = delta_log
  )

make_log2fc_heatmap_for_ora <- function(ct, ct_res,
                                        max_genes = MAX_GENES_PER_CT,
                                        base = HEAT_BASE_SIZE,
                                        return_plot = FALSE) {
  
  if (!nrow(ct_res)) return(NULL)
  
  df <- purrr::map_dfr(seq_len(nrow(ct_res)), function(i) {
    genes_i <- ct_res$overlapGenes[[i]]
    if (!length(genes_i)) return(tibble())
    
    tibble(
      CellType   = ct,
      pathway    = ct_res$pathway[i],
      path_index = ct_res$path_index[i],
      Gene       = genes_i,
      sig        = ct_res$sig[i]
    )
  })
  
  if (!nrow(df)) return(NULL)
  
  keep_genes <- df %>%
    count(Gene, sort = TRUE, name = "n_path") %>%
    left_join(
      df %>%
        group_by(Gene) %>%
        summarise(max_sig = max(sig, na.rm = TRUE), .groups = "drop"),
      by = "Gene"
    ) %>%
    arrange(desc(n_path), desc(max_sig), Gene) %>%
    slice_head(n = max_genes) %>%
    pull(Gene)
  
  plot_df <- df %>%
    filter(Gene %in% keep_genes) %>%
    left_join(
      logfc_long %>% filter(CellType == ct),
      by = c("CellType", "Gene")
    )
  
  if (!nrow(plot_df)) return(NULL)
  
  gene_levels <- plot_df %>%
    group_by(Gene) %>%
    summarise(
      n_path = n_distinct(pathway),
      max_sig = max(sig, na.rm = TRUE),
      mean_abs_diff = mean(abs(log2FC_chronic - log2FC_nonchronic), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_path), desc(mean_abs_diff), desc(max_sig), Gene) %>%
    pull(Gene)
  
  heat_df <- plot_df %>%
    select(path_index, Gene, log2FC_chronic, log2FC_nonchronic) %>%
    pivot_longer(
      cols = c(log2FC_nonchronic, log2FC_chronic),
      names_to = "Group",
      values_to = "log2FC"
    ) %>%
    mutate(
      Group_short = recode(
        Group,
        log2FC_nonchronic = "NC",
        log2FC_chronic = "C"
      ),
      Gene = factor(Gene, levels = gene_levels),
      y = if_else(
        Group_short == "NC",
        path_index + 0.22,
        path_index - 0.22
      )
    ) %>%
    arrange(path_index, desc(Group_short))
  
  if (!nrow(heat_df)) return(NULL)
  
  group_breaks <- c(rbind(ct_res$path_index + 0.22, ct_res$path_index - 0.22))
  group_labels <- rep(c("NC", "C"), times = nrow(ct_res))
  
  # Dynamic symmetric color scale per cell type
  heat_lim <- quantile(abs(heat_df$log2FC), 0.95, na.rm = TRUE)
  
  if (!is.finite(heat_lim) || heat_lim == 0) {
    heat_lim <- 1
  }
  
  p <- ggplot(heat_df, aes(x = Gene, y = y, fill = log2FC)) +
    geom_tile(color = "grey90", linewidth = 0.25, height = 0.30) +
    geom_hline(
      yintercept = seq(1.5, nrow(ct_res) - 0.5, by = 1),
      color = "grey85",
      linewidth = 0.35
    ) +
    scale_fill_gradient2(
    low = "steelblue2",
    mid = "grey95",
    high = "firebrick2",
    midpoint = 0,
    limits = c(-heat_lim, heat_lim),
    oob = scales::squish,
    na.value = "white",
    name = "log2FC A→6M"
  ) +
    scale_y_continuous(
      breaks = group_breaks,
      labels = group_labels,
      limits = c(0.5, nrow(ct_res) + 0.5),
      expand = c(0, 0)
    ) +
    labs(x = "Overlapping genes", y = NULL) +
    theme_classic(base_size = base) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.6),
      axis.text.x = element_text(
        size = GENE_TEXT_SIZE,
        angle = GENE_LABEL_ANGLE,
        hjust = 1,
        vjust = 1,
        margin = margin(t = 3)
      ),
      axis.text.y = element_text(size = GROUP_TEXT_SIZE),
      legend.title = element_text(size = LEGEND_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      plot.margin = margin(6, 4, 6, 0)
    )
  
  if (return_plot) return(p)
  invisible(p)
}

make_dotplot_for_shared_ora <- function(ct_res,
                                        max_genes = MAX_GENES_PER_CT,
                                        base = HEAT_BASE_SIZE,
                                        return_plot = FALSE) {
  
  if (!nrow(ct_res)) return(NULL)
  
  dot_df <- purrr::map_dfr(seq_len(nrow(ct_res)), function(i) {
    genes_i <- ct_res$overlapGenes[[i]]
    if (!length(genes_i)) return(tibble())
    tibble(
      path_index = ct_res$path_index[i],
      Gene       = genes_i,
      effect     = ct_res$effect[i],
      sig        = ct_res$sig[i]
    )
  })
  
  if (!nrow(dot_df)) return(NULL)
  
  keep_genes <- dot_df |>
    count(Gene, sort = TRUE, name = "n_path") |>
    left_join(
      dot_df |>
        group_by(Gene) |>
        summarise(max_sig = max(sig, na.rm = TRUE), .groups = "drop"),
      by = "Gene"
    ) |>
    arrange(desc(n_path), desc(max_sig)) |>
    slice_head(n = max_genes) |>
    pull(Gene)
  
  dot_df <- dot_df |>
    filter(Gene %in% keep_genes) |>
    mutate(Gene = factor(Gene, levels = keep_genes))
  
  if (!nrow(dot_df)) return(NULL)
  
  sig_max <- max(dot_df$sig, na.rm = TRUE)
  eff_max <- max(abs(dot_df$effect), na.rm = TRUE)
  if (!is.finite(sig_max) || sig_max == 0) sig_max <- 1
  if (!is.finite(eff_max) || eff_max == 0) eff_max <- 1
  
  p <- ggplot(dot_df, aes(x = Gene, y = path_index)) +
    geom_point(aes(size = sig, color = effect), alpha = 0.9) +
    scale_size_continuous(
      name = "-log10\nadj. p",
      range = c(2.5, 9),
      limits = c(0, sig_max)
    ) +
    scale_color_gradient2(
      low = "steelblue2",
      mid = "grey88",
      high = "firebrick2",
      midpoint = 0,
      limits = c(-eff_max, eff_max),
      name = "ORA\neffect"
    ) +
    scale_y_continuous(
      breaks = ct_res$path_index,
      labels = rep("", nrow(ct_res)),
      limits = c(0.5, nrow(ct_res) + 0.5),
      expand = c(0, 0)
    ) +
    labs(x = "Overlapping gene", y = NULL) +
    theme_minimal(base_size = base) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size = GENE_TEXT_SIZE,
        angle = GENE_LABEL_ANGLE,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = LEGEND_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      plot.margin = margin(6, 4, 6, 0)
    )
  
  if (return_plot) return(p)
  invisible(p)
}

compose_panel <- function(ct, ora_table, out_subdir, subtitle_prefix = "Reactome ORA") {
  
  n_total_pathways <- ora_table %>%
    filter(CellType == ct, !is.na(padj_u)) %>%
    nrow()
  
  ct_res <- ora_table %>%
    filter(CellType == ct, !is.na(padj_u)) %>%
    arrange(padj_u) %>%
    slice_head(n = ifelse(ct == shared_group_label, collapse_bar_n, bar_n)) %>%
    mutate(
      label = str_wrap(pretty_pathway(pathway), width = WRAP_WIDTH_PATH),
      sig = safe_sig(padj_u),
      path_index = row_number()
    ) %>%
    filter(is.finite(sig), is.finite(effect))
  
  if (!nrow(ct_res)) return(NULL)
  
  sig_max <- max(ct_res$sig, na.rm = TRUE)
  if (!is.finite(sig_max) || sig_max == 0) sig_max <- 1
  
  panel_height <- max(5, 0.55 * nrow(ct_res) + 3)
  
  p_bar <- ggplot(ct_res, aes(x = sig, y = path_index, fill = sig)) +
    geom_col(width = 0.55, orientation = "y") +
    scale_y_continuous(
      breaks = ct_res$path_index,
      labels = ct_res$label,
      limits = c(0.5, nrow(ct_res) + 0.5),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = c(0, sig_max * 1.1),
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    scale_fill_gradient(
      low = "steelblue2",
      high = "firebrick2",
      name = "-log10\nadj. p",
      limits = c(0, 15),
      oob = scales::squish
    ) +
    labs(x = "-log10 adjusted p-value", y = NULL) +
    theme_classic(base_size = BAR_BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.6),
      axis.text.y = element_text(size = PATHWAY_TEXT_SIZE, lineheight = 1.25, hjust = 1),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      legend.title = element_text(size = LEGEND_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      plot.margin = margin(6, 10, 6, 10)
    ) +
    coord_cartesian(clip = "off")
  
  if (ct == shared_group_label) {
    p_right <- make_dotplot_for_shared_ora(
      ct_res = ct_res,
      max_genes = MAX_GENES_PER_CT,
      base = HEAT_BASE_SIZE,
      return_plot = TRUE
    )
  } else {
    p_right <- make_log2fc_heatmap_for_ora(
      ct = ct,
      ct_res = ct_res,
      max_genes = MAX_GENES_PER_CT,
      base = HEAT_BASE_SIZE,
      return_plot = TRUE
    )
  }
  
  if (is.null(p_right)) {
    p_right <- patchwork::plot_spacer()
  }
  
  panel <- p_bar + p_right +
    patchwork::plot_layout(widths = c(1.45, 3.35)) +
    patchwork::plot_annotation(
      title = pretty_ct_ora(ct),
      subtitle = paste0(
        subtitle_prefix,
        "  |  n = ", n_total_pathways, " pathways",
        "  |  |Δlog2FC| ≥ ", delta_threshold,
        "  |  minSize = ", minSize_ora,
        "  |  maxSize = ", maxSize_ora,
        "  |  FDR ≤ ", ora_fdr_cutoff
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = TITLE_SIZE_PT, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30")
      )
    )
  
  base_name <- file.path(output_dir, out_subdir, paste0("panel_", safe(ct)))
  
  ggsave(
    paste0(base_name, ".png"),
    panel,
    width = PANEL_W_IN,
    height = panel_height,
    dpi = 300
  )
  
  tryCatch(
    ggsave(
      paste0(base_name, ".pdf"),
      panel,
      width = PANEL_W_IN,
      height = panel_height,
      device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
    ),
    error = function(e) message("PDF save failed for ", ct, ": ", e$message)
  )
  
  invisible(panel)
}


unlink(Sys.glob(file.path(output_dir, "reactome_panels_ora_heatmap", "*.png")))
unlink(Sys.glob(file.path(output_dir, "reactome_panels_ora_heatmap", "*.pdf")))
unlink(Sys.glob(file.path(output_dir, "reactome_panels_ora_heatmap_collapsed", "*.png")))
unlink(Sys.glob(file.path(output_dir, "reactome_panels_ora_heatmap_collapsed", "*.pdf")))
# ============================================================
# RENDER ORA PANELS
# ============================================================

cts_to_plot <- c(ora_targets, shared_group_label)
cts_to_plot <- unique(cts_to_plot)
cts_to_plot <- intersect(cts_to_plot, unique(ORA_results$CellType))

if (length(cts_to_plot) > 0 && nrow(ORA_results) > 0) {
  message("\nPlotting ", length(cts_to_plot), " standard ORA panel(s)...")
  
  panels_standard <- purrr::map(
    cts_to_plot,
    ~ compose_panel(
      ct = .x,
      ora_table = ORA_results,
      out_subdir = "reactome_panels_ora_heatmap",
      subtitle_prefix = "Reactome ORA"
    )
  )
  panels_standard <- purrr::compact(panels_standard)
  
  if (length(panels_standard)) {
    pdf(
      file.path(output_dir, "reactome_panels_ora_heatmap", "all_panels_ora_heatmap.pdf"),
      width = PANEL_W_IN,
      height = PANEL_H_IN,
      useDingbats = FALSE
    )
    purrr::walk(panels_standard, print)
    dev.off()
  }
}

if (collapse_enabled && length(cts_to_plot) > 0 && nrow(ORA_results_collapsed) > 0) {
  message("\nPlotting ", length(cts_to_plot), " collapsed ORA panel(s)...")
  
  panels_collapsed <- purrr::map(
    cts_to_plot,
    ~ compose_panel(
      ct = .x,
      ora_table = ORA_results_collapsed,
      out_subdir = "reactome_panels_ora_heatmap_collapsed",
      subtitle_prefix = "Reactome ORA (collapsed)"
    )
  )
  panels_collapsed <- purrr::compact(panels_collapsed)
  
  if (length(panels_collapsed)) {
    pdf(
      file.path(output_dir, "reactome_panels_ora_heatmap_collapsed", "all_panels_ora_heatmap_collapsed.pdf"),
      width = PANEL_W_IN,
      height = PANEL_H_IN,
      useDingbats = FALSE
    )
    purrr::walk(panels_collapsed, print)
    dev.off()
  }
}



# ============================================================
# COMBINE MAIN FIGURE ON A4 PAGE
# ============================================================

# Extract plot objects needed for main figure
p_scatter_nonclassical <- scatter_plots[["Non_classical_monocytes"]]
p_scatter_pdc <- scatter_plots[["Plasmacytoid_dendritic_cells"]]

# Store ORA panels as objects
p_ora_nonclassical <- compose_panel(
  ct = "Non_classical_monocytes",
  ora_table = ORA_results_collapsed,
  out_subdir = "reactome_panels_ora_heatmap_collapsed",
  subtitle_prefix = "Reactome ORA (collapsed)"
)

p_ora_pdc <- compose_panel(
  ct = "Plasmacytoid_dendritic_cells",
  ora_table = ORA_results_collapsed,
  out_subdir = "reactome_panels_ora_heatmap_collapsed",
  subtitle_prefix = "Reactome ORA (collapsed)"
)


library(patchwork)

# A4 landscape dimensions in inches
A4_W <- 11.69
A4_H <- 8.27

# Example object names — replace with your actual plot objects
# p_bar_summary
# p_scatter_nonclassical
# p_scatter_pdc
# p_ora_nonclassical
# p_traj_nonclassical
# p_ora_pdc
# p_traj_pdc

main_figure <- 
  (
    p_bar_summary + 
      (p_scatter_nonclassical / p_scatter_pdc)
  ) +
  (
    p_ora_nonclassical / 
      p_ora_pdc
  ) +
  (
    p_traj_nonclassical / 
      p_traj_pdc
  ) +
  plot_layout(
    widths = c(1.05, 1.35, 1.05)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave(
  filename = file.path(output_dir, "main_figure_A4_landscape.png"),
  plot = main_figure,
  width = A4_W,
  height = A4_H,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "main_figure_A4_landscape.pdf"),
  plot = main_figure,
  width = A4_W,
  height = A4_H,
  bg = "white",
  limitsize = FALSE
)
# ============================================================
# SECOND VERSION OF MAIN FIGURE:
# ORA significance bubbles + curated trajectory plots only
# No heatmap
# Does NOT remove the original figure
# ============================================================

library(patchwork)

A4_W <- 11.69
A4_H <- 8.27

p_scatter_nonclassical <- scatter_plots[["Non_classical_monocytes"]]
p_scatter_pdc <- scatter_plots[["Plasmacytoid_dendritic_cells"]]

out_v2_dir <- file.path(output_dir, "main_figure_v2_bubble_trajectories_only")

dir.create(
  out_v2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ----------------------------
# Curated genes to show in Figure 6C/D
# ----------------------------

nonclassical_selected_genes <- c(
  "SH3RF3",
  "SMAD3",
  "CDK14",
  "HIPK2",
  "FNIP1",
  "HLA-DPA1",
  "HLA-DPB1",
  "FCGR3A",
  "RHOA",
  "CFL1"
)

pdc_selected_genes <- c(
  "HIVEP3",
  "IFI44L",
  "KMT2C",
  "NFKB2",
  "CDYL",
  "TCF4",
  "LILRA4",
  "GZMB",
  "HLA-DPA1",
  "CST3"
)

intermediate_selected_genes <- c(
  "IFI44L",
  "HIPK2",
  "PARP14",
  "ANPEP",
  "LIMK2",
  "HLA-DPA1",
  "HLA-DPB1",
  "FCGR3A",
  "KLF2",
  "SLC25A5"
)

# ----------------------------
# Shorten long Reactome pathway labels manually
# ----------------------------

shorten_pathway_label <- function(x) {
  x <- pretty_pathway(x)
  
  x <- gsub(
    "ANTIGEN PRESENTATION FOLDING ASSEMBLY AND PEPTIDE LOADING OF CLASS I MHC",
    "CLASS I MHC\nANTIGEN PRESENTATION FOLDING AND PEPTIDE LOADING",
    x
  )
  
  x <- gsub(
    "SMAD2 SMAD3 SMAD4 HETEROTRIMER REGULATES TRANSCRIPTION",
    "SMAD2/3/4\nREGULATES TRANSCRIPTION",
    x
  )
  
  x <- gsub(
    "IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL",
    "IMMUNOREGULATORY\nLYMPHOID-NON-LYMPHOID\nINTERACTIONS",
    x
  )
  
  x <- gsub(
    "PLATELET ACTIVATION SIGNALING AND AGGREGATION",
    "PLATELET ACTIVATION,\nSIGNALING AND AGGREGATION",
    x
  )
  
  x <- gsub(
    "FCGAMMA RECEPTOR FCGR DEPENDENT PHAGOCYTOSIS",
    "FCgammaR-DEPENDENT\nPHAGOCYTOSIS",
    x
  )
  
  x <- gsub(
    "RHO GTPASES ACTIVATE IQGAPS",
    "RHO GTPASES\nACTIVATE IQGAPS",
    x
  )
  
  x
}

# ----------------------------
# ORA bubble plot
# x-axis, size and fill = -log10 adjusted p-value
# number inside dot = number of overlapping genes
# ----------------------------

make_ora_bubble_plot_compact <- function(ct, ora_table, top_n = 10, base_size = 8.5) {
  
  ct_res <- ora_table %>%
    filter(CellType == ct, !is.na(padj_u)) %>%
    arrange(padj_u) %>%
    slice_head(n = top_n) %>%
    mutate(
      Pathway = vapply(pathway, shorten_pathway_label, character(1)),
      Pathway = str_wrap(Pathway, width =36),
      sig = -log10(padj_u),
      Pathway = factor(Pathway, levels = rev(Pathway))
    ) %>%
    filter(is.finite(sig))
  
  if (!nrow(ct_res)) return(plot_spacer())
  
  sig_max <- max(ct_res$sig, na.rm = TRUE)
  if (!is.finite(sig_max) || sig_max == 0) sig_max <- 1
  
  ggplot(ct_res, aes(x = sig, y = Pathway)) +
    geom_point(
      aes(size = sig, fill = sig),
      shape = 21,
      colour = "grey25",
      stroke = 0.35,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = overlap),
      size = 2.5,
      colour = "black"
    ) +
    scale_size_continuous(
      name = "-log10\nadj. p",
      range = c(2.8, 7.2),
      limits = c(0, sig_max)
    ) +
    scale_fill_gradient(
      low = "steelblue2",
      high = "firebrick2",
      guide = "none",
      limits = c(0, sig_max),
      oob = scales::squish
    ) +
    scale_x_continuous(
      limits = c(0, sig_max * 1.05),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "-log10 adjusted p-value",
      y = NULL,
      caption = "Numbers inside bubbles = overlapping genes"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      aspect.ratio = 1.22,
      axis.text.y = element_text(
        size = 7.2,
        lineheight = 0.92,
        margin = margin(r = 5)
      ),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.35, "cm"),
      plot.caption = element_text(size = 7.5, hjust = 0.5, margin = margin(t = 4)),
      plot.margin = margin(4, 14, 4, 4)
    ) +
    coord_cartesian(clip = "off")
}


# ----------------------------
# Curated-gene trajectory plot
# Facet labels include only gene name + Δlog(C-NC)
# ----------------------------

make_curated_trajectory_plot_compact <- function(ct, curated_genes, ncol = 5, base_size = 7.2) {
  
  curated_genes <- toupper(curated_genes)
  
  if (!length(curated_genes)) return(plot_spacer())
  
  f_C_main      <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))
  f_NC_main     <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))
  f_C_sig_in_NC <- file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))
  f_NC_sig_in_C <- file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv"))
  
  C_main_expr   <- read_deg_expr_full(f_C_main, "C")
  NC_main_expr  <- read_deg_expr_full(f_NC_main, "NC")
  C_other_expr  <- read_deg_expr_full(f_C_sig_in_NC, "C")
  NC_other_expr <- read_deg_expr_full(f_NC_sig_in_C, "NC")
  
  expr_merged <- full_join(
    C_main_expr %>%
      rename(
        C_A_expr_main = C_A_expr,
        C_M6_expr_main = C_M6_expr,
        C_log2FC_main = C_log2FC
      ),
    C_other_expr %>%
      rename(
        C_A_expr_other = C_A_expr,
        C_M6_expr_other = C_M6_expr,
        C_log2FC_other = C_log2FC
      ),
    by = "gene"
  ) %>%
    full_join(
      NC_main_expr %>%
        rename(
          NC_A_expr_main = NC_A_expr,
          NC_M6_expr_main = NC_M6_expr,
          NC_log2FC_main = NC_log2FC
        ),
      by = "gene"
    ) %>%
    full_join(
      NC_other_expr %>%
        rename(
          NC_A_expr_other = NC_A_expr,
          NC_M6_expr_other = NC_M6_expr,
          NC_log2FC_other = NC_log2FC
        ),
      by = "gene"
    ) %>%
    mutate(
      Gene = toupper(gene),
      C_A_expr   = coalesce(C_A_expr_main, C_A_expr_other),
      C_M6_expr  = coalesce(C_M6_expr_main, C_M6_expr_other),
      NC_A_expr  = coalesce(NC_A_expr_main, NC_A_expr_other),
      NC_M6_expr = coalesce(NC_M6_expr_main, NC_M6_expr_other),
      C_log2FC   = coalesce(C_log2FC_main, C_log2FC_other),
      NC_log2FC  = coalesce(NC_log2FC_main, NC_log2FC_other),
      delta_log  = C_log2FC - NC_log2FC
    ) %>%
    filter(Gene %in% curated_genes) %>%
    mutate(
      gene_order = match(Gene, curated_genes)
    ) %>%
    arrange(gene_order)
  
  missing_genes <- setdiff(curated_genes, expr_merged$Gene)
  if (length(missing_genes) > 0) {
    warning(
      "These curated genes were not found for ", ct, ": ",
      paste(missing_genes, collapse = ", ")
    )
  }
  
  if (!nrow(expr_merged)) return(plot_spacer())
  
  label_levels <- expr_merged %>%
    mutate(
      Gene_label = paste0(
        Gene,
        "\n\u0394log(C-NC) = ", round(delta_log, 2)
      )
    ) %>%
    arrange(gene_order) %>%
    pull(Gene_label)
  
  expr_merged <- expr_merged %>%
    mutate(
      Gene_label = paste0(
        Gene,
        "\n\u0394log(C-NC) = ", round(delta_log, 2)
      ),
      Gene_label = factor(Gene_label, levels = label_levels)
    )
  
  plot_df <- expr_merged %>%
    pivot_longer(
      cols = c(C_A_expr, C_M6_expr, NC_A_expr, NC_M6_expr),
      names_to = "group_time",
      values_to = "avg_expr"
    ) %>%
    mutate(
      Outcome = case_when(
        str_starts(group_time, "C_")  ~ "Chronic",
        str_starts(group_time, "NC_") ~ "Non-chronic"
      ),
      Time = case_when(
        str_detect(group_time, "_A_expr$")  ~ "Acute",
        str_detect(group_time, "_M6_expr$") ~ "6M"
      ),
      Time = factor(Time, levels = c("Acute", "6M")),
      Outcome = factor(Outcome, levels = c("Non-chronic", "Chronic"))
    ) %>%
    filter(!is.na(avg_expr))
  
  ggplot(plot_df, aes(x = Time, y = avg_expr, group = Outcome, colour = Outcome)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3) +
    facet_wrap(~ Gene_label, scales = "free_y", ncol = ncol) +
    scale_colour_manual(values = c(
      "Non-chronic" = "#4575b4",
      "Chronic" = "#d73027"
    )) +
    labs(
      x = NULL,
      y = "Average expression",
      colour = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 8, lineheight = 0.95),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ----------------------------
# Compose one ORA + trajectory panel
# Includes pathway-analysis specifics in subtitle
# Also explains bubble numbers and Δlog meaning
# ----------------------------

compose_panel_v2 <- function(ct, ora_table, curated_genes, top_n_pathways = 10) {
  
  curated_genes <- toupper(curated_genes)
  
  ct_all <- ora_table %>%
    filter(CellType == ct, !is.na(padj_u))
  
  n_total_pathways <- nrow(ct_all)
  
  p_bubble <- make_ora_bubble_plot_compact(
    ct = ct,
    ora_table = ora_table,
    top_n = top_n_pathways,
    base_size = 8.5
  )
  
  p_traj <- make_curated_trajectory_plot_compact(
    ct = ct,
    curated_genes = curated_genes,
    ncol = 5,
    base_size = 6.6
  )
  
  panel <- p_bubble + p_traj +
    plot_layout(widths = c(1.10, 2.05)) +
    plot_annotation(
      title = pretty_celltype(ct),
      subtitle = paste0(
        "Reactome ORA (collapsed) | n = ", n_total_pathways,
        " | |\u0394log2FC| >= ", delta_threshold,
        " | minSize = 5 | maxSize = 500 | FDR <= 0.05"
            ),
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 7.6, colour = "grey30", lineheight = 1.05)
      )
    )
  
  return(panel)
}

# ----------------------------
# Build C/D panels for main figure v2
# ----------------------------

message("Building v2 non-classical monocyte panel...")

p_ora_nonclassical_v2 <- compose_panel_v2(
  ct = "Non_classical_monocytes",
  ora_table = ORA_results_collapsed,
  curated_genes = nonclassical_selected_genes,
  top_n_pathways = 10
)

message("Building v2 pDC panel...")

p_ora_pdc_v2 <- compose_panel_v2(
  ct = "Plasmacytoid_dendritic_cells",
  ora_table = ORA_results_collapsed,
  curated_genes = pdc_selected_genes,
  top_n_pathways = 10
)

message("Building v2 intermediate monocyte supplementary panel...")

p_ora_intermediate_v2 <- compose_panel_v2(
  ct = "Intermediate_monocytes",
  ora_table = ORA_results_collapsed,
  curated_genes = intermediate_selected_genes,
  top_n_pathways = 10
)

# ----------------------------
# Cleaner scatter panels for B
# ----------------------------

p_scatter_nonclassical_v2 <- p_scatter_nonclassical +
  labs(title = "Non-classical monocytes") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

p_scatter_pdc_v2 <- p_scatter_pdc +
  labs(title = "pDCs") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

# ----------------------------
# Main figure v2
# ----------------------------

main_figure_v2 <- 
  (
    p_bar_summary +
      (p_scatter_nonclassical_v2 / p_scatter_pdc_v2)
  ) /
  (
    p_ora_nonclassical_v2 /
      p_ora_pdc_v2
  ) +
  plot_layout(
    heights = c(0.8, 1.55)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

# ----------------------------
# Save outputs
# ----------------------------

message("Saving main figure v2...")

out_main_png <- file.path(
  out_v2_dir,
  "main_figure_A4_landscape_v2_bubble_trajectories_only.png"
)

out_main_pdf <- file.path(
  out_v2_dir,
  "main_figure_A4_landscape_v2_bubble_trajectories_only.pdf"
)

ggsave(
  filename = out_main_png,
  plot = main_figure_v2,
  width = A4_W,
  height = A4_H,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Saved: ", out_main_png)

ggsave(
  filename = out_main_pdf,
  plot = main_figure_v2,
  width = A4_W,
  height = A4_H,
  bg = "white",
  limitsize = FALSE
)

message("Saved: ", out_main_pdf)

message("Saving separate C/D panels...")

out_C <- file.path(out_v2_dir, "panel_C_nonclassical_v2.png")
out_D <- file.path(out_v2_dir, "panel_D_pDC_v2.png")

ggsave(
  filename = out_C,
  plot = p_ora_nonclassical_v2,
  width = 13,
  height = 4.2,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Saved: ", out_C)

ggsave(
  filename = out_D,
  plot = p_ora_pdc_v2,
  width = 13,
  height = 4.2,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Saved: ", out_D)

out_IM <- file.path(out_v2_dir, "supp_panel_intermediate_monocytes_v2.png")

ggsave(
  filename = out_IM,
  plot = p_ora_intermediate_v2,
  width = 13,
  height = 4.2,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Saved: ", out_IM)
# ============================================================
# FINAL DIAGNOSTICS
# ============================================================

cat("=== 1. FILE COMPLETENESS CHECK ===\n")
file_check <- map_dfr(celltypes, function(ct) {
  tibble(
    cell_type     = ct,
    C_main        = file.exists(file.path(input_dir, paste0("C.A_vs_6M.", ct, ".tsv"))),
    NC_main       = file.exists(file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".tsv"))),
    C_sig_in_NC   = file.exists(file.path(input_dir, paste0("C.A_vs_6M.", ct, ".sig_in_NC.tsv"))),
    NC_sig_in_C   = file.exists(file.path(input_dir, paste0("NC.A_vs_6M.", ct, ".sig_in_C.tsv")))
  )
}) %>%
  mutate(complete = C_main & NC_main & C_sig_in_NC & NC_sig_in_C)

print(file_check)
cat("\nIncomplete cell types:\n")
print(file_check %>% filter(!complete))

cat("\n=== 2. COLUMN NAME / DESIRED_ORDER CHECK ===\n")
remaining_cols_check <- setdiff(
  colnames(shared_heatmap_df)[colnames(shared_heatmap_df) != "Gene"],
  desired_order
)
if (length(remaining_cols_check) == 0) {
  cat("OK: all heatmap columns are covered by desired_order\n")
} else {
  cat("WARNING: these columns are NOT in desired_order and will be appended unordered:\n")
  print(remaining_cols_check)
}

unmatched_desired <- desired_order[!desired_order %in% colnames(shared_heatmap_df)]
if (length(unmatched_desired) == 0) {
  cat("OK: all desired_order entries match actual column names\n")
} else {
  cat("WARNING: these desired_order entries do not match any column (typo?):\n")
  print(unmatched_desired)
}

cat("\n=== 3. PRETTY_CELLTYPE CONSISTENCY CHECK ===\n")
pretty_names <- pretty_celltype(celltypes)
unprettied <- pretty_names[pretty_names == celltypes]
if (length(unprettied) == 0) {
  cat("OK: all cell types were transformed by pretty_celltype()\n")
} else {
  cat("WARNING: these cell types were NOT changed by pretty_celltype() - may be missing a rule:\n")
  print(unprettied)
}

cat("\n=== 4. DELTA LOG SANITY CHECK ===\n")
cat("Rows in combined: ", nrow(combined), "\n")
cat("NA delta_log: ", sum(is.na(combined$delta_log)), "\n")
cat("Cell types in combined: ", n_distinct(combined$cell_type), "\n")
cat("Expected cell types: ", length(celltypes), "\n")
missing_from_combined <- setdiff(celltypes, unique(combined$cell_type))
if (length(missing_from_combined) == 0) {
  cat("OK: all cell types present in combined\n")
} else {
  cat("WARNING: these cell types are missing from combined (were skipped in main loop):\n")
  print(missing_from_combined)
}

cat("\n=== 5. ORA INPUT CHECK ===\n")
cat("Universe size: ", length(tested_gene_universe), "\n")
cat("Pathways after universe filtering: ", length(gs_react), "\n")
cat("Genes >= ", shared_group_cutoff, " groups: ", length(genes_ge4), "\n")
cat("ORA targets: ", paste(ora_targets, collapse = ", "), "\n")
missing_ora <- setdiff(ora_targets, unique(long_ora$CellType))
if (length(missing_ora) == 0) {
  cat("OK: all ORA targets present in long_ora\n")
} else {
  cat("WARNING: these ORA targets are missing from long_ora:\n")
  print(missing_ora)
}

cat("\n=== 6. ORA RESULTS CHECK ===\n")
if (nrow(ORA_results) == 0) {
  cat("WARNING: ORA_results is empty - no pathways passed filters\n")
} else {
  cat("ORA results summary:\n")
  print(
    ORA_results %>%
      group_by(CellType) %>%
      summarise(
        n_pathways = n(),
        min_padj   = min(padj_u, na.rm = TRUE),
        max_padj   = max(padj_u, na.rm = TRUE),
        .groups    = "drop"
      )
  )
}

cat("\n=== 7. TRAJECTORY FORCE INCLUDE CHECK ===\n")
missing_forced <- setdiff(trajectory_force_include, celltypes)
if (length(missing_forced) == 0) {
  cat("OK: all forced trajectory cell types exist in detected celltypes\n")
} else {
  cat("WARNING: these forced cell types were not detected (wrong name?):\n")
  print(missing_forced)
}

cat("\n=== ALL DIAGNOSTICS COMPLETE ===\n")

message("ALL DONE.")