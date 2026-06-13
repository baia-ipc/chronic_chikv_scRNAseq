# =========================================================
# Reactome pathway overview — Non-classical monocytes, Acute phase
# Single cell type, single file
# Themes: 3 enriched in chronic, 10 depleted in chronic
# Vertical version: pathways on y-axis
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(cowplot)
  library(ggnewscale)
})

# ---------------------------------------------------------------
# User settings
# ---------------------------------------------------------------
data_dir <- "F:/chikv-results/reanalysis/pathways"

files <- list(
  `Non-clas. mono` = "TabV06_reactome_gsea_1_25_Non_classical_monocytes.A.C_vs_NC.Non_classical_monocytes.GSEA.Reactome.tsv"
)

out_dir <- file.path(data_dir, "reactome_overview_nonclassical_mono_acute")
dir.create(out_dir, showWarnings = FALSE)

padj_cutoff <- 0.05
cell_order  <- c("Non-clas. mono")

# ---------------------------------------------------------------
# Theme definitions
# ---------------------------------------------------------------
theme_map <- tribble(
  ~theme_key,                ~display_label,                     ~regex,
  
  "NOTCH signaling",         "NOTCH signaling",
  "notch4 intracellular|runx3 regulates notch",
  
  "Chromatin organization",  "Chromatin organization",
  "^chromatin organization$|^chromatin modifying enzymes$|epigenetic regulation of gene expression|mll3 and mll4|wdr5-containing histone",
  
  "RHO GTPase cycles",      "RHO/RAC GTPase cycling",
  "^rac1 gtpase cycle$|^rho gtpase cycle$",
  
  "RHO effectors",          "RHO effectors (ROCK/CIT)",
  "rho gtpases activate rocks|rho gtpases activate rhotekin|rho gtpases activate cit",
  
  "Interferon signaling",    "Interferon signaling",
  "^interferon signaling$|interferon alpha/beta signaling|^interferon gamma signaling$",
  
  "Mitochondrial respiration","Mitochondrial respiration",
  "respiratory electron transport|aerobic respiration|formation of atp by chemiosmotic coupling|complex i biogenesis|complex iv assembly|cristae formation|complex iii assembly",
  
  "Mitochondrial translation","Mitochondrial translation",
  "^mitochondrial translation$|mitochondrial translation elongation|mitochondrial translation initiation|mitochondrial translation termination",
  
  "Translation elongation",  "Translation / elongation",
  "^eukaryotic translation elongation$|^translation elongation$|^translation termination$|^translation$|srp-dependent cotranslational",
  
  "Antigen presentation",    "Antigen presentation",
  "^antigen processing-cross presentation$|^er-phagosome pathway$|cross-presentation of soluble exogenous antigens|antigen processing: ub, atp-independent|mhc class ii antigen presentation|class i mhc mediated",
  
  "NF-kB signaling",         "NF-kB signaling",
  "nf-kappab in b cells|tnfr2 non-canonical|nik.*nf-kb|fceri mediated nf-kb|dectin-1 mediated noncanonical nf-kb",
  
  "Cytokine signaling",      "Cytokine signaling",
  "^cytokine signaling in immune system$|^signaling by interleukins$",
  
  "Cell cycle",              "Cell cycle",
  "^cell cycle$|^cell cycle, mitotic$|^g1/s transition$|^s phase$|^mitotic g2|^g2/m transition$|^m phase$|apc/c|^dna replication$|dna replication pre-initiation",
  
  "DNA damage TP53",         "DNA damage / TP53",
  "tp53 regulates|transcriptional regulation by tp53|p53-dependent|p53-independent|stabilization of p53|g1/s dna damage|g2 checkpoint",
  
  "Biological oxidations",   "Biological oxidations / Phase II",
  "^biological oxidations$|^phase ii - conjugation of compounds$|^glutathione conjugation$"
)

assign_theme_key <- function(x) {
  x2 <- tolower(x)
  hit <- theme_map$theme_key[
    map_lgl(theme_map$regex, ~ str_detect(x2, regex(.x, ignore_case = TRUE)))
  ]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ---------------------------------------------------------------
# Read file
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
# Row/column order
# ---------------------------------------------------------------
theme_keys_topdown <- c(
  "NOTCH signaling",
  "Chromatin organization",
  "RHO GTPase cycles",        # ← up (was RAC1 RHO GTPase)
  "Interferon signaling",
  "Mitochondrial respiration",
  "Mitochondrial translation",
  "Translation elongation",
  "Antigen presentation",
  "NF-kB signaling",
  "RHO effectors",            # ← down (new)
  "Cytokine signaling",
  "Cell cycle",
  "DNA damage TP53",
  "Biological oxidations"
)

n_up   <- 3
n_down <- 11

direction_lookup <- setNames(
  c(rep("up", n_up), rep("down", n_down)),
  theme_keys_topdown
)

label_lookup <- setNames(theme_map$display_label, theme_map$theme_key)

# ---------------------------------------------------------------
# Build plot grid
# ---------------------------------------------------------------
plot_grid_df <- expand.grid(
  theme_key  = theme_keys_topdown,
  cell_group = cell_order,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(best_per_theme, by = c("theme_key", "cell_group")) %>%
  mutate(
    present     = !is.na(NES),
    significant = if_else(is.na(significant), FALSE, significant),
    absNES      = if_else(is.na(NES), 0, abs(NES)),
    log10padj   = if_else(is.na(padj) | padj == 0, 0, -log10(padj)),
    direction   = direction_lookup[theme_key]
  )

# ---------------------------------------------------------------
# Theme audit
# ---------------------------------------------------------------
gsea_audit <- gsea_all %>%
  mutate(
    theme_key   = map_chr(description, assign_theme_key),
    significant = !is.na(padj) & padj <= padj_cutoff
  ) %>%
  filter(!is.na(theme_key)) %>%
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
  left_join(best_flags, by = c("theme", "cell", "pathway"), suffix = c("", ".best")) %>%
  mutate(is_best = if_else(is.na(is_best.best), is_best, is_best.best)) %>%
  select(-is_best.best)

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
    theme_key  = theme_keys_topdown,
    cell_group = cell_order,
    fill = list(
      n_matched = 0,
      n_significant = 0,
      best_pathway = NA,
      best_NES = NA,
      best_padj = NA
    )
  )

write_csv(audit_detail,  file.path(out_dir, "theme_audit_detail_mono_acute.csv"))
write_csv(audit_summary, file.path(out_dir, "theme_audit_summary_mono_acute.csv"))
message("Audit tables saved.")

# ---------------------------------------------------------------
# Add n_significant to plot_grid
# ---------------------------------------------------------------
n_sig_lookup <- audit_summary %>%
  select(theme_key, cell_group, n_significant) %>%
  mutate(
    theme_key  = as.character(theme_key),
    cell_group = as.character(cell_group)
  )

plot_grid_df <- plot_grid_df %>%
  mutate(
    theme_key  = as.character(theme_key),
    cell_group = as.character(cell_group)
  ) %>%
  select(-any_of("n_significant")) %>%
  left_join(n_sig_lookup, by = c("theme_key", "cell_group")) %>%
  mutate(
    n_significant = tidyr::replace_na(n_significant, 0)
  )

# ---------------------------------------------------------------
# Plot (vertical: pathways on y-axis)
# ---------------------------------------------------------------
plot_grid_df <- plot_grid_df %>%
  mutate(
    theme_key  = factor(theme_key, levels = rev(theme_keys_topdown)),
    cell_group = factor(cell_group, levels = cell_order)
  )

bg_data <- tibble(
  theme_key = factor(theme_keys_topdown, levels = rev(theme_keys_topdown)),
  direction = direction_lookup[theme_keys_topdown],
  x_mid     = 1,
  width     = 0.96
)

# center labels on actual colored blocks
label_pos <- bg_data %>%
  mutate(y_num = as.numeric(theme_key)) %>%
  group_by(direction) %>%
  summarise(y = mean(y_num), .groups = "drop")

y_up   <- label_pos$y[label_pos$direction == "up"]
y_down <- label_pos$y[label_pos$direction == "down"]

nes_lim <- ceiling(max(abs(plot_grid_df$NES), na.rm = TRUE) * 4) / 4

p_main <- ggplot() +
  geom_tile(
    data = bg_data,
    aes(x = x_mid, y = theme_key, width = width, height = 0.96, fill = direction),
    alpha = 0.12, color = NA, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("up" = "#d73027", "down" = "#2b6cb0"),
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_point(
    data = plot_grid_df %>% filter(!present | !significant),
    aes(x = cell_group, y = theme_key),
    shape = 21, size = 6, fill = "white", color = "grey72", stroke = 0.8
  ) +
  
  geom_point(
    data = plot_grid_df %>% filter(present, significant),
    aes(x = cell_group, y = theme_key, size = log10padj, fill = NES),
    shape = 21, color = "grey25", stroke = 0.2, alpha = 0.95
  ) +
  
  scale_fill_gradientn(
    colours = c("#2b6cb0", "#6baed6", "#f7f7f7", "#fc8d59", "#d73027"),
    values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
    limits  = c(-3, 3),
    oob     = scales::squish,
    name    = "NES"
  ) +
  
  scale_size_continuous(
    range = c(6, 13),
    name  = expression(-log[10](p[adj]))
  ) +
  
  scale_size_continuous(range = c(6, 13), name = "|NES|") +
  
  geom_text(
    data = plot_grid_df %>% filter(present, significant),
    aes(x = cell_group, y = theme_key, label = n_significant),
    size = 3.2, color = "white", fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  scale_y_discrete(labels = label_lookup[theme_keys_topdown]) +
  scale_x_discrete(
    position = "top",
    expand = expansion(mult = c(0.55, 0.05))
  ) + 
  
  coord_cartesian(clip = "off") +
  
  labs(
    x       = NULL,
    y       = NULL,
    title   = NULL,
    caption = "Dot size = -log10(adj. p-value). Dot colour = NES. \nNumber = significant pathways per theme."
  )+
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(face = "bold", size = 15, color = "black"),
    axis.text.y        = element_text(size = 13, color = "black"),
    legend.position    = "none",
    plot.caption       = element_text(size = 11, color = "grey40",
                                      hjust = 0, margin = margin(t = 10)),
    plot.margin        = margin(20, 10, 50, 80)
  ) +
  annotate(
    "text", x = 0.48, y = y_up,
    label = "\u2191 enriched in chronic",
    angle = 90, hjust = 0.5, vjust = 0.5,
    color = "#d73027", fontface = "bold", size = 4.5
  ) +
  annotate(
    "text", x = 0.48, y = y_down,
    label = "\u2193 depleted in chronic",
    angle = 90, hjust = 0.5, vjust = 0.5,
    color = "#2b6cb0", fontface = "bold", size = 4.5
  )

# ---------------------------------------------------------------
# Extract legend
# ---------------------------------------------------------------
p_for_legend <- p_main +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black")
  )

leg <- cowplot::get_legend(p_for_legend)

# ---------------------------------------------------------------
# Combine main plot + legend
# ---------------------------------------------------------------
final <- cowplot::plot_grid(
  p_main, leg,
  nrow = 1,
  rel_widths = c(0.78, 0.22)
)

ggsave(file.path(out_dir, "Reactome_overview_nonclassical_mono_acute.pdf"),
       final, width = 7.5, height = 7.5, bg = "white")
ggsave(file.path(out_dir, "Reactome_overview_nonclassical_mono_acute.png"),
       final, width = 7.5, height = 7.5, dpi = 300, bg = "white")

message("Done. Output in: ", out_dir)
final