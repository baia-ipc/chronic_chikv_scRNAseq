# =========================================================
# 6-month grouped DEG heatmap with ComplexHeatmap
# - brighter colors
# - clearer log2FC legend
# - module legend below log2FC legend
# - genes ordered within each module by strongest cell group
# - genes with |log2FC| < 0.25 masked
# - only one set of legends
# - no module names printed on the left
# - narrower columns
# - includes newly added genes
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -----------------------------
# User settings
# -----------------------------
file_path  <- "F://chikv-results//reanalysis//DEGs.xlsx"
sheet_name <- "6 months "
out_dir    <- "figure_6m_grouped_heatmap"
dir.create(out_dir, showWarnings = FALSE)

wanted_cells <- c("PBMC", "Naive CD4", "pDC", "NK", "Th17", "Th1/17")
min_abs_log2fc <- 0.25

# -----------------------------
# Helpers
# -----------------------------
read_deg_sheet <- function(file_path, sheet_name) {
  raw <- read_excel(file_path, sheet = sheet_name, col_names = FALSE)
  
  celltypes <- as.character(unlist(raw[1, ]))
  celltypes <- str_trim(celltypes)
  dat <- raw[-c(1, 2), ]
  pair_starts <- seq(1, ncol(dat), by = 2)
  
  map_dfr(pair_starts, function(i) {
    if ((i + 1) > ncol(dat)) return(NULL)
    ct <- celltypes[i]
    if (is.na(ct) || ct == "") return(NULL)
    
    tibble(
      gene = dat[[i]],
      log2FC = dat[[i + 1]],
      celltype = ct
    ) %>%
      mutate(
        gene = as.character(gene),
        gene = str_trim(gene),
        celltype = as.character(celltype),
        celltype = str_trim(celltype),
        log2FC = as.character(log2FC),
        log2FC = str_replace_all(log2FC, ",", "."),
        log2FC = suppressWarnings(as.numeric(log2FC))
      ) %>%
      filter(!is.na(gene), gene != "", !is.na(log2FC))
  })
}

shorten_celltypes <- function(x) {
  recode(
    x,
    "Naive CD4 T" = "Naive CD4",
    "pDCs" = "pDC",
    "NK cells" = "NK",
    "Th17" = "Th17",
    "Th1/Th17" = "Th1/17",
    "total cells" = "PBMC",
    .default = x
  )
}

# -----------------------------
# Read + filter
# -----------------------------
deg <- read_deg_sheet(file_path, sheet_name) %>%
  mutate(celltype = shorten_celltypes(celltype)) %>%
  filter(celltype %in% wanted_cells)

# Keep genes with >= 0.25 in at least one selected cell type
deg6 <- deg %>%
  group_by(gene) %>%
  filter(any(abs(log2FC) >= min_abs_log2fc)) %>%
  ungroup()

# -----------------------------
# Build matrix
# -----------------------------
mat <- deg6 %>%
  distinct(gene, celltype, .keep_all = TRUE) %>%
  pivot_wider(names_from = celltype, values_from = log2FC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

wanted_cells_present <- wanted_cells[wanted_cells %in% colnames(mat)]
mat <- mat[, wanted_cells_present, drop = FALSE]

# Mask weak individual values
mat_masked <- mat
mat_masked[abs(mat_masked) < min_abs_log2fc] <- NA

# -----------------------------
# Curated modules
# -----------------------------
gene_module_tbl <- tibble(
  gene = rownames(mat_masked),
  module = case_when(
    gene %in% c("GZMB", "KIR3DL2", "KLRF1", "PRF1", "SIGLEC7", "SPON2", "NCR1") ~ "Cytotoxic / NK effector",
    
    gene %in% c("ANXA1", "ARID5A", "BACH2", "CD200", "CITED2", "IL12RB1", "IRF4", "LGALS1", "TCF7", "TOX2", "PTGDS", "FCRL3") ~ "Immune regulation / differentiation",
    
    gene %in% c("ABCA1", "ABCG1", "GAPDH", "GSTP1", "MT1F", "MTHFD1L", "SIK1B", "SLC25A5") ~ "Metabolism / mitochondria",
    
    gene %in% c("ABR", "ARPC1B", "CAMK1D", "CCR7", "CD82", "LSP1", "MICAL3", "PIK3CA", "PRKCH", "PTPRE", "RAB35", "RALA", "SPTBN1", "TTC7A", "RHOC", "MOB1B") ~ "Signaling / trafficking",
    
    gene %in% c("AUTS2", "DACH1", "DDX21", "GTF2IRD2", "H3F3A", "KLF7", "MAML3", "TET3", "ZNF704") ~ "Transcription / chromatin",
    
    gene %in% c("BBC3", "BIRC2", "N4BP1", "DDI2") ~ "Stress / apoptosis",
    
    gene %in% c("ADAM23", "CLIC3", "GRN", "GLT8D1", "HDDC3", "ISM1", "TMEM123", "TMEM41B", "YIPF5", "VSIG1") ~ "Membrane / secretory biology",
    
    gene %in% c("RADX") ~ "DNA damage / genome stability",
    
    gene %in% c("AC025164.1", "LINC00299", "SNHG29", "MIR181A2HG", "AP005019.1") ~ "Non-coding / lncRNA",
    
    TRUE ~ "Unassigned"
  )
) %>%
  filter(module != "Unassigned")

mat_masked <- mat_masked[rownames(mat_masked) %in% gene_module_tbl$gene, , drop = FALSE]

# -----------------------------
# Order genes within module by strongest cell group
# -----------------------------
module_levels <- c(
  "Cytotoxic / NK effector",
  "Immune regulation / differentiation",
  "Metabolism / mitochondria",
  "Signaling / trafficking",
  "Transcription / chromatin",
  "Stress / apoptosis",
  "Membrane / secretory biology",
  "DNA damage / genome stability",
  "Non-coding / lncRNA"
)

gene_module_tbl <- gene_module_tbl %>%
  mutate(module = factor(module, levels = module_levels))

gene_order_tbl <- as.data.frame(mat_masked) %>%
  rownames_to_column("gene") %>%
  rowwise() %>%
  mutate(
    strongest_cell = {
      vals <- c_across(all_of(wanted_cells_present))
      if (all(is.na(vals))) NA_character_ else wanted_cells_present[which.max(abs(vals))]
    },
    strongest_value = {
      vals <- c_across(all_of(wanted_cells_present))
      if (all(is.na(vals))) NA_real_ else max(abs(vals), na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  mutate(strongest_cell = factor(strongest_cell, levels = wanted_cells_present))

gene_module_tbl <- gene_module_tbl %>%
  left_join(gene_order_tbl %>% select(gene, strongest_cell, strongest_value), by = "gene")

ordered_genes <- gene_module_tbl %>%
  arrange(module, strongest_cell, desc(strongest_value), gene) %>%
  pull(gene)

mat_ord <- mat_masked[ordered_genes, , drop = FALSE]

row_split <- gene_module_tbl %>%
  slice(match(rownames(mat_ord), gene)) %>%
  pull(module)

row_split <- factor(row_split, levels = module_levels)

# -----------------------------
# Colors
# -----------------------------
max_abs <- max(abs(mat_ord), na.rm = TRUE)
lim <- max(1.5, ceiling(max_abs * 4) / 4)

# Brighter blue-white-red
col_fun <- colorRamp2(
  c(-lim, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, lim),
  c("#0033FF", "#3366FF", "#7FAAFF", "#DCE8FF", "white",
    "#FFD9D9", "#FF8C8C", "#FF3333", "#CC0000")
)

module_cols <- c(
  "Cytotoxic / NK effector" = "#E33022",
  "Immune regulation / differentiation" = "#9B4FB3",
  "Metabolism / mitochondria" = "#1FA77E",
  "Signaling / trafficking" = "#FF8700",
  "Transcription / chromatin" = "#53B541",
  "Stress / apoptosis" = "#B55E24",
  "Membrane / secretory biology" = "#63C3AD",
  "DNA damage / genome stability" = "#F6D32D",
  "Non-coding / lncRNA" = "#A9A9A9"
)

# -----------------------------
# Left module bar
# -----------------------------
module_factor <- factor(
  row_split,
  levels = module_levels
)

left_anno <- rowAnnotation(
  module = module_factor,
  col = list(module = module_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(5, "mm")
)
# -----------------------------
# Legends
# -----------------------------
lgd_fc <- Legend(
  title = "log2FC (C vs NC)",
  col_fun = col_fun,
  at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
  labels = c("↓ Chronic", "-1.5", "-1.0", "-0.5", "0", "0.5", "1.0", "1.5", "↑ Chronic"),
  direction = "vertical",
  legend_height = unit(6, "cm"),
  title_gp = gpar(fontsize = 11, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

lgd_module <- Legend(
  title = "module",
  at = names(module_cols),
  legend_gp = gpar(fill = module_cols, col = NA),
  ncol = 1,
  title_gp = gpar(fontsize = 11, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

combined_lgd <- packLegend(
  lgd_fc,
  lgd_module,
  direction = "vertical",
  gap = unit(6, "mm")
)

# -----------------------------
# Heatmap
# -----------------------------
ht <- Heatmap(
  mat_ord,
  name = "log2FC",
  col = col_fun,
  na_col = "grey95",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = row_split,
  row_gap = unit(2, "mm"),
  column_names_rot = 45,
  
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  rect_gp = gpar(col = "#D0D0D0", lwd = 0.6),
  
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  
  left_annotation = left_anno,
  
  # narrower columns
  width = unit(4.2, "cm"),
  
  # do not print module names on left
  row_title = NULL,
  
  # remove built-in legends
  show_heatmap_legend = FALSE
)

# -----------------------------
# Draw and save
# -----------------------------
pdf(file.path(out_dir, "6m_grouped_heatmap_complex.pdf"), width = 7.8, height = 11, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_list = list(combined_lgd),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)
dev.off()

png(file.path(out_dir, "6m_grouped_heatmap_complex.png"), width = 7.8, height = 11, units = "in", res = 300)
draw(
  ht,
  heatmap_legend_list = list(combined_lgd),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)
dev.off()

# -----------------------------
# Save outputs
# -----------------------------
write.csv(
  cbind(gene = rownames(mat_ord), as.data.frame(mat_ord)),
  file.path(out_dir, "6m_grouped_heatmap_matrix.csv"),
  row.names = FALSE
)

write.csv(
  gene_module_tbl,
  file.path(out_dir, "6m_grouped_heatmap_modules.csv"),
  row.names = FALSE
)