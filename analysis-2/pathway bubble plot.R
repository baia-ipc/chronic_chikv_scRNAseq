# =========================================================
# Reactome pathway overview — v5
# Fill now shows -log10(adjusted p-value), like your
# non-classical monocyte script
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(ggnewscale)
})

data_dir <- "F:/chikv-results/reanalysis/pathways"

files <- list(
  pDC         = "TabV06_reactome_gsea_2_20_Plasmacytoid_dendritic_cells.6M.C_vs_NC.Plasmacytoid_dendritic_cells.GSEA.Reactome.tsv",
  NK          = "TabV06_reactome_gsea_2_2_Natural_killer_cells.6M.C_vs_NC.Natural_killer_cells.GSEA.Reactome.tsv",
  `Naive CD4` = "TabV06_reactome_gsea_2_3_Naive_CD4_T_cells.6M.C_vs_NC.Naive_CD4_T_cells.GSEA.Reactome.tsv",
  Th17        = "TabV06_reactome_gsea_2_9_Th17_cells.6M.C_vs_NC.Th17_cells.GSEA.Reactome.tsv",
  `Th1/17`    = "TabV06_reactome_gsea_2_12_Th1_Th17_cells.6M.C_vs_NC.Th1_Th17_cells.GSEA.Reactome.tsv"
)

out_dir <- file.path(data_dir, "reactome_overview_v5")
dir.create(out_dir, showWarnings = FALSE)

padj_cutoff <- 0.05
cell_order  <- c("pDC", "NK", "Naive CD4", "Th17", "Th1/17")

# ---------------------------------------------------------------
# Theme definitions — plain ASCII keys only
# ---------------------------------------------------------------
theme_map <- tribble(
  ~theme_key,                   ~display_label,                              ~regex,
  
  "Interferon signaling",  "Interferon signaling",
  "^interferon signaling$|interferon gamma signaling|ddx58/ifih1-mediated induction of interferon",
  
  "Mitochondrial translation",   "Mitochondrial translation",
  "^mitochondrial translation$|mitochondrial translation elongation|mitochondrial translation initiation|mitochondrial translation termination",
  
  "Mitochondrial respiration",   "Mitochondrial respiration",
  "^respiratory electron transport$|aerobic respiration and respiratory electron transport|^formation of atp by chemiosmotic coupling$|^complex i biogenesis$|^complex iv assembly$|^cristae formation$",
  
  "Translation elongation",      "Translation / elongation",
  "^eukaryotic translation elongation$|^translation elongation$|^translation termination$|^translation initiation$|^translation$",
  
  "Chaperonin CCT folding",      "Chaperonin / CCT-TriC folding",
  "cooperation of prefoldin|prefoldin mediated transfer|formation of tubulin folding intermediates|post-chaperonin tubulin folding|folding of actin by cct",
  
  "Platelet Ca2 degranulation",  "Platelet / Ca\u00b2\u207a degranulation",
  "^platelet degranulation$|^response to elevated platelet cytosolic ca2",
  
  "Detoxification ROS",          "Detoxification / ROS",
  "^detoxification of reactive oxygen species$",
  
  "Biological oxidations",       "Biological oxidations / Phase II",
  "^biological oxidations$|^phase ii - conjugation of compounds$|^glutathione conjugation$|^cytosolic sulfonation of small molecules$",
  
  "Antigen presentation",        "Antigen presentation",
  "^antigen processing-cross presentation$|^er-phagosome pathway$|^cross-presentation of soluble exogenous antigens|^antigen processing: ub, atp-independent",
  
  "TGF-beta SMAD signaling",     "TGF-\u03b2 / SMAD signaling",
  "^signaling by tgf-beta receptor complex$|^signaling by tgfb family members$|^signaling by tgfbr3$|^tgfbr3 expression$|^downregulation of smad2|^transcriptional activity of smad2",
  
  "NR1H2 NR1H3 signaling",       "NR1H2 / NR1H3 signaling",
  "nr1h2 and nr1h3-mediated signaling|nr1h3 & nr1h2 regulate gene expression",
  
  "Chromatin organization",      "Chromatin organization",
  "^chromatin organization$|^chromatin modifying enzymes$|^epigenetic regulation of gene expression$",
  
  "SUMOylation",                 "SUMOylation",
  "^sumoylation$|sumoylation of dna|sumoylation of rna|sumoylation of sumoylation|sumoylation of ubiquitin",
  
  "Pre-NOTCH transcription",     "Pre-NOTCH transcription",
  "pre-notch transcription",
  
  "Nuclear pore disassembly",    "Nuclear pore complex disassembly",
  "^nuclear pore complex.*disassembly$|postmitotic nuclear pore",
  
  "Circadian clock",             "Circadian clock",
  "^circadian clock$|bmal1.*activates circadian|expression of bmal|phosphorylated bmal1|cry:per:kinase",
  
  "Cytokine signaling",          "Cytokine signaling",
  "^cytokine signaling in immune system$|^signaling by interleukins$"
)

assign_theme_key <- function(x) {
  x2  <- tolower(x)
  hit <- theme_map$theme_key[
    map_lgl(theme_map$regex, ~ str_detect(x2, regex(.x, ignore_case = TRUE)))
  ]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ---------------------------------------------------------------
# Read files
# ---------------------------------------------------------------
read_gsea_tsv <- function(path, cell_group) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    transmute(
      description = Description,
      NES         = as.numeric(NES),
      padj        = as.numeric(`p.adjust`),
      cell_group  = cell_group
    ) %>%
    filter(!is.na(NES))
}

gsea_all <- imap_dfr(files, ~ read_gsea_tsv(file.path(data_dir, .x), .y))

gsea_themed <- gsea_all %>%
  mutate(
    theme_key   = map_chr(description, assign_theme_key),
    significant = !is.na(padj) & padj <= padj_cutoff
  ) %>%
  filter(!is.na(theme_key))

best_per_theme <- gsea_themed %>%
  group_by(cell_group, theme_key) %>%
  slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup()

message("Themes matched: ", paste(sort(unique(best_per_theme$theme_key)), collapse = ", "))

# ---------------------------------------------------------------
# Row order — exactly as you want TOP to BOTTOM
# ---------------------------------------------------------------
theme_keys_topdown <- c(
  "Mitochondrial respiration",
  "Mitochondrial translation",
  "Translation elongation",
  "Chaperonin CCT folding",
  "Platelet Ca2 degranulation",
  "Detoxification ROS",
  "Biological oxidations",
  "Antigen presentation",
  "NR1H2 NR1H3 signaling",
  "Chromatin organization",
  "SUMOylation",
  "Pre-NOTCH transcription",
  "TGF-beta SMAD signaling",
  "Nuclear pore disassembly",
  "Circadian clock",
  "Interferon signaling",
  "Cytokine signaling"
)

n_up   <- 8
n_down <- 9

direction_lookup <- setNames(
  c(rep("up", n_up), rep("down", n_down)),
  theme_keys_topdown
)

label_lookup <- setNames(theme_map$display_label, theme_map$theme_key)

# ---------------------------------------------------------------
# Build plot grid
# ---------------------------------------------------------------
theme_keys_levels <- rev(theme_keys_topdown)

plot_grid <- expand.grid(
  theme_key  = theme_keys_topdown,
  cell_group = cell_order,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(best_per_theme, by = c("theme_key", "cell_group")) %>%
  mutate(
    present      = !is.na(NES),
    significant  = if_else(is.na(significant), FALSE, significant),
    absNES       = if_else(is.na(NES), 0, abs(NES)),
    log10padj    = if_else(is.na(padj) | padj <= 0, 0, -log10(padj)),
    direction    = direction_lookup[theme_key],
    theme_key    = factor(theme_key,  levels = theme_keys_levels),
    cell_group   = factor(cell_group, levels = cell_order)
  )

message("NA theme_key after factor: ", sum(is.na(plot_grid$theme_key)))

# ---------------------------------------------------------------
# Background data
# ---------------------------------------------------------------
bg_data <- tibble(
  theme_key = factor(theme_keys_topdown, levels = theme_keys_levels),
  direction = direction_lookup[theme_keys_topdown],
  x_mid     = mean(seq_along(cell_order)),
  width     = length(cell_order) + 0.96
)

# Side annotation y positions
up_keys   <- theme_keys_topdown[1:n_up]
down_keys <- theme_keys_topdown[(n_up + 1):length(theme_keys_topdown)]

y_up   <- mean(match(up_keys,   theme_keys_levels))
y_down <- mean(match(down_keys, theme_keys_levels))

# max for p-adjust legend
padj_fill_max <- max(plot_grid$log10padj, na.rm = TRUE)
if (!is.finite(padj_fill_max) || padj_fill_max == 0) padj_fill_max <- 1

# ---------------------------------------------------------------
# Theme audit
# ---------------------------------------------------------------
gsea_audit <- gsea_all %>%
  mutate(
    theme_key   = map_chr(description, assign_theme_key),
    significant = !is.na(padj) & padj <= padj_cutoff
  ) %>%
  filter(!is.na(theme_key)) %>%
  mutate(
    cell_group = factor(cell_group, levels = cell_order),
    theme_key  = factor(theme_key,  levels = theme_keys_topdown)
  ) %>%
  arrange(theme_key, cell_group, desc(abs(NES)))

audit_detail <- gsea_audit %>%
  transmute(
    theme       = theme_key,
    cell        = cell_group,
    pathway     = description,
    NES         = round(NES, 3),
    padj        = round(padj, 4),
    significant = significant,
    is_best     = FALSE
  )

best_flags <- gsea_audit %>%
  group_by(theme_key, cell_group) %>%
  slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(theme = theme_key, cell = cell_group, pathway = description, is_best = TRUE)

audit_detail <- audit_detail %>%
  rows_update(best_flags, by = c("theme", "cell", "pathway"))

audit_summary <- gsea_audit %>%
  group_by(theme_key, cell_group) %>%
  summarise(
    n_matched     = n(),
    n_significant = sum(significant),
    best_pathway  = description[which.max(abs(NES))],
    best_NES      = round(NES[which.max(abs(NES))], 3),
    best_padj     = round(padj[which.max(abs(NES))], 4),
    .groups = "drop"
  ) %>%
  complete(
    theme_key  = factor(theme_keys_topdown, levels = theme_keys_topdown),
    cell_group = factor(cell_order, levels = cell_order),
    fill = list(
      n_matched = 0,
      n_significant = 0,
      best_pathway = NA,
      best_NES = NA,
      best_padj = NA
    )
  )

# ---------------------------------------------------------------
# Console audit
# ---------------------------------------------------------------
walk(theme_keys_topdown, function(tk) {
  cat(sprintf("\n=== %s ===\n", tk))
  sub <- audit_detail %>% filter(theme == tk)
  if (nrow(sub) == 0) {
    cat("  (no matches in any cell type)\n")
    return()
  }
  walk(cell_order, function(cl) {
    csub <- sub %>% filter(cell == cl)
    if (nrow(csub) == 0) {
      cat(sprintf("  %-12s  (no match)\n", cl))
    } else {
      walk(seq_len(nrow(csub)), function(i) {
        best_flag <- if (csub$is_best[i]) " [PLOTTED]" else ""
        sig_flag  <- if (csub$significant[i]) "*" else " "
        cat(sprintf("  %-12s  NES=%+.3f  padj=%.4f%s  %s%s\n",
                    cl, csub$NES[i], csub$padj[i],
                    sig_flag, csub$pathway[i], best_flag))
      })
    }
  })
})

cat("\n\n=== SUMMARY: n significant pathways per theme x cell ===\n")
audit_summary %>%
  select(theme_key, cell_group, n_matched, n_significant) %>%
  pivot_wider(
    names_from  = cell_group,
    values_from = c(n_matched, n_significant),
    names_glue  = "{cell_group}_{.value}"
  ) %>%
  print(n = Inf, width = Inf)

write_csv(audit_detail,  file.path(out_dir, "theme_audit_detail.csv"))
write_csv(audit_summary, file.path(out_dir, "theme_audit_summary.csv"))
message("Audit tables saved to: ", out_dir)

# ---------------------------------------------------------------
# Add n_significant to plot_grid
# ---------------------------------------------------------------
n_sig_lookup <- audit_summary %>%
  select(theme_key, cell_group, n_significant) %>%
  mutate(
    theme_key  = factor(theme_key,  levels = theme_keys_levels),
    cell_group = factor(cell_group, levels = cell_order)
  )

plot_grid <- plot_grid %>%
  left_join(n_sig_lookup, by = c("theme_key", "cell_group")) %>%
  mutate(
    n_significant = replace_na(n_significant, 0)
  )

# ---------------------------------------------------------------
# Plot
# ---------------------------------------------------------------

nes_abs_max <- max(abs(plot_grid$NES[plot_grid$present & plot_grid$significant]),
                   na.rm = TRUE)
nes_abs_max <- ceiling(nes_abs_max * 10) / 10   # round up to 1 dp

p <- ggplot() +
  
  geom_tile(
    data = bg_data,
    aes(x = x_mid, y = theme_key, width = width, height = 0.96, fill = direction),
    alpha = 0.12, color = NA, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("up" = "#d73027", "down" = "#2b6cb0"),
    guide  = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # non-significant / absent: hollow circles
  geom_point(
    data  = plot_grid %>% filter(!present | !significant),
    aes(x = cell_group, y = theme_key),
    shape = 21, size = 4.5, fill = "white", color = "grey72", stroke = 0.8
  ) +
  
  # significant dots — SIZE = -log10(padj), FILL = NES
  geom_point(
    data  = plot_grid %>% filter(present, significant),
    aes(x = cell_group, y = theme_key,
        size = log10padj,        # <-- was absNES
        fill = NES),             # <-- was log10padj
    shape = 21, color = "grey25", stroke = 0.2, alpha = 0.95
  ) +
  
  # Diverging fill for NES (blue = negative / depleted, red = positive / enriched)
  scale_fill_gradientn(
    colours = c("#2b6cb0", "#6baed6", "#f7f7f7", "#fc8d59", "#d73027"),
    values  = scales::rescale(c(-2.8, -1.6, 0, 1.6, 2.8), 
                              from = c(-3, 3)),
    limits  = c(-3, 3),
    oob     = scales::squish,
    name    = "NES"
  ) +
  
  # Size now encodes -log10(padj)
  scale_size_continuous(
    range = c(6, 13),
    name  = expression(-log[10](p[adj]))
  ) +
  
  geom_text(
    data = plot_grid %>% filter(present, significant),
    aes(x = cell_group, y = theme_key, label = n_significant),
    size = 3, color = "grey20", fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  scale_y_discrete(labels = label_lookup[theme_keys_levels]) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  
  labs(
    x       = NULL,
    y       = NULL,
    title   = NULL,
    caption = "Dot size = -log10(adj. p-value). Dot colour = NES. Number = significant pathways per theme."
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x.top    = element_text(face = "bold", size = 16, color = "black"),
    axis.text.y        = element_text(size = 15, color = "black"),
    plot.title         = element_text(face = "bold", size = 15, color = "black"),
    legend.title       = element_text(face = "bold", color = "black"),
    legend.text        = element_text(color = "black"),
    plot.caption       = element_text(
      size = 12, color = "black",
      hjust = 0, margin = margin(t = 10)
    ),
    plot.margin        = margin(20, 100, 20, 100)
  ) +
  annotate(
    "text", x = 0.1, y = y_up,
    label = "\u2191 enriched in chronic", angle = 90,
    hjust = 0.5, vjust = 0.5, color = "#d73027", fontface = "bold", size = 5
  ) +
  annotate(
    "text", x = 0.1, y = y_down,
    label = "\u2193 depleted in chronic", angle = 90,
    hjust = 0.5, vjust = 0.5, color = "#2b6cb0", fontface = "bold", size = 5
  )

ggsave(file.path(out_dir, "Reactome_overview_v6.pdf"),
       p, width = 13, height = 12.5, bg = "white")
ggsave(file.path(out_dir, "Reactome_overview_v6.png"),
       p, width = 13, height = 12.5, dpi = 300, bg = "white")

message("Done. Output in: ", out_dir)
p