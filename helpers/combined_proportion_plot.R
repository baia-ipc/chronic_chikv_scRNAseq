library(cowplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

proportion_barplot <- function(
    seurat_obj,
    celltype_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    point_color_col = NULL,
    point_color_title = NULL,
    point_color_palette = NULL,
    log_scale = F
) {
  if (is.null(point_color_col)) {
    # default to coloring by condition_col
    point_color_col <- condition_col
    point_color_title <- ifelse(is.null(point_color_title),
                                condition_title, point_color_title)
    plot_data <- seurat_obj@meta.data %>%
      group_by(sample, .data[[condition_col]], .data[[celltype_col]]) %>%
      summarize(cell_count = n(), .groups = 'drop') %>%
      group_by(sample, .data[[condition_col]]) %>%
      mutate(fraction = cell_count / sum(cell_count) * 100) %>%
      ungroup()
  } else {
    # if user provided a different column, we still need a title
    if (is.null(point_color_title)) point_color_title <- point_color_col
    plot_data <- seurat_obj@meta.data %>%
      group_by(sample, .data[[condition_col]], .data[[celltype_col]],
               .data[[point_color_col]]) %>%
      summarize(cell_count = n(), .groups = 'drop') %>%
      group_by(sample, .data[[condition_col]]) %>%
      mutate(fraction = cell_count / sum(cell_count) * 100) %>%
      ungroup()
  }

  custom_color <- condition_cols
  names(custom_color) <- levels(seurat_obj@meta.data[[condition_col]])
  custom_fill <- scales::alpha(condition_cols, 0.5)
  names(custom_fill) <- levels(seurat_obj@meta.data[[condition_col]])

  # optionally define palette for point_color_col
  if (is.null(point_color_palette)) {
    # if point_color_col == condition_col, use original palette
    if (point_color_col == condition_col) {
      point_color_palette <- custom_color
    } else {
      # create a default discrete palette for any other coloring column
      unique_vals <- unique(seurat_obj@meta.data[[point_color_col]])
      default_cols <-
        RColorBrewer::brewer.pal(min(length(unique_vals), 9), "Set1")
      point_color_palette <-
        setNames(default_cols[seq_along(unique_vals)], unique_vals)
    }
  }

  p <- ggplot(plot_data, aes_string(x = celltype_col, y = "fraction")) +
    geom_boxplot(aes_string(fill = condition_col), alpha = 0.5,
                 outlier.shape = NA, position = position_dodge(width = 0.9)) +
    geom_point(aes_string(color = point_color_col, group = condition_col),
               position = position_dodge(width = 0.9), size = 1.6) +
    scale_fill_manual(values = custom_fill, guide = 'none') +
    scale_color_manual(values = point_color_palette) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.border = element_rect(size = 1, fill = NA),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                 size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x.bottom = element_text(size = 10, color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 14, color = "black"),
      strip.text.x = element_text(angle = 0, color = "black"),
      strip.text.y = element_text(angle = 0, color = "black"),
    ) +
    labs(
      x = NULL,
      y = "Proportion of cells (%)",
      fill = condition_title,
      color = point_color_title
    )

  if (log_scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}

grouped_proportion_barplot <- function(
    seurat_obj,
    celltype_col,
    cell_group_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    panel_bg_cols = NULL,
    point_color_col = NULL,
    point_color_title = NULL,
    point_color_palette = NULL,
    log_scale = F
) {
  if (is.null(point_color_col)) {
    # default to coloring by condition_col
    point_color_col <- condition_col
    point_color_title <- ifelse(is.null(point_color_title), condition_title, point_color_title)
    plot_data <- seurat_obj@meta.data %>%
      group_by(sample, .data[[condition_col]], .data[[celltype_col]], .data[[cell_group_col]]) %>%
      summarize(cell_count = n(), .groups = 'drop') %>%
      group_by(sample, .data[[condition_col]]) %>%
      mutate(fraction = cell_count / sum(cell_count) * 100) %>%
      ungroup()
  } else {
    # if user provided a different column, we still need a title
    if (is.null(point_color_title)) point_color_title <- point_color_col
    plot_data <- seurat_obj@meta.data %>%
      group_by(sample, .data[[condition_col]], .data[[celltype_col]], .data[[cell_group_col]], .data[[point_color_col]]) %>%
      summarize(cell_count = n(), .groups = 'drop') %>%
      group_by(sample, .data[[condition_col]]) %>%
      mutate(fraction = cell_count / sum(cell_count) * 100) %>%
      ungroup()
  }

  bg_data <- unique(plot_data[, c(cell_group_col, celltype_col)])
  bg_data[[condition_col]] <- "dummy"
  bg_data$fraction <- 1

  # color/fill for condition
  custom_color <- condition_cols
  names(custom_color) <- levels(seurat_obj@meta.data[[condition_col]])
  custom_fill <- scales::alpha(condition_cols, 0.5)
  names(custom_fill) <- levels(seurat_obj@meta.data[[condition_col]])

  # optionally define palette for point_color_col
  if (is.null(point_color_palette)) {
    # if point_color_col == condition_col, use original palette
    if (point_color_col == condition_col) {
      point_color_palette <- custom_color
    } else {
      # create a default discrete palette for any other coloring column
      unique_vals <- unique(seurat_obj@meta.data[[point_color_col]])
      default_cols <-
        RColorBrewer::brewer.pal(min(length(unique_vals), 9), "Set1")
      point_color_palette <-
        setNames(default_cols[seq_along(unique_vals)], unique_vals)
    }
  }

  if (is.null(panel_bg_cols)) {
    unique_groups <- length(unique(seurat_obj@meta.data[[cell_group_col]]))
    panel_bg_cols <- RColorBrewer::brewer.pal(min(unique_groups, 9), "Pastel1")
    names(panel_bg_cols) <- unique(seurat_obj@meta.data[[cell_group_col]])
  } else {
    panel_bg_cols <- scales::alpha(panel_bg_cols, 0.2)
  }

  p <- ggplot(plot_data, aes_string(x = celltype_col, y = "fraction")) +
    geom_rect(
      data = bg_data,
      aes_string(fill = cell_group_col),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2
    ) +
    scale_fill_manual(values = c(custom_fill, panel_bg_cols), guide = 'none') +
    geom_boxplot(aes_string(fill = condition_col), alpha = 0.5,
                 outlier.shape = NA, position = position_dodge(width = 0.9)) +
    geom_point(
      aes_string(color = point_color_col, group = condition_col),
      position = position_dodge(width = 0.9),
      size = 1.6
    ) +
    scale_color_manual(values = point_color_palette) +
    facet_grid(as.formula(paste("~", cell_group_col)),
               scales = "free_x", space = "free") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.border = element_rect(size = 1, fill = NA),
      axis.text.y = element_text(size = 11.5, color = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                 size = 11.5, color = "black"),
      strip.text = element_text(size = 14.5, color = "black"),
      strip.text.y = element_text(angle = 0, color = "black")
    ) +
    labs(
      x = NULL,
      y = "Proportion of cells (%)",
      fill = condition_title,
      color = point_color_title
    )

  if (log_scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}

permutation_test_plot <- function(
    sc_utils_obj,
    celltype_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    point_color_by_celltype = F,
    point_color_title = NULL,
    point_color_palette = NULL
) {
  meta_data <- as.data.frame(sc_utils_obj@meta_data)
  permutation_data <- as.data.frame(sc_utils_obj@results$permutation)

  plot_data <- permutation_data %>%
    left_join(meta_data, by = c("clusters" = celltype_col)) %>%
    mutate(significance =
             ifelse(FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
                    paste("FDR <", FDR_threshold, "&", 
                          expression("|Log2(FD)| >="),
                          round(log2FD_threshold, 2)), "n.s."),
           significance =
             factor(significance,
                    levels = c(paste("FDR <", FDR_threshold, "&",
                                     expression("|Log2(FD)| >="),
                                     round(log2FD_threshold, 2)), "n.s.")))

  plot_data$clusters <- factor(plot_data$clusters,
                               levels = levels(meta_data[[celltype_col]]))
  if (!point_color_by_celltype) {
    point_color_col = "significance"
    point_color_palette = c("salmon", "grey")
    point_color_title = "Significance"
  } else {
    point_color_col <- "clusters"
    if (is.null(point_color_palette)) {
      unique_vals <- unique(plot_data$clusters)
      default_cols <-
        RColorBrewer::brewer.pal(min(length(unique_vals), 9), "Set1")
      point_color_palette <-
        setNames(default_cols[seq_along(unique_vals)], unique_vals)
    }
    if (is.null(point_color_title)) {
      point_color_title <- "Celltype"
    }
  }

  p <- ggplot(plot_data, aes_string(x = "clusters", y = "obs_log2FD")) +
    geom_pointrange(aes_string(ymin = "boot_CI_2.5", ymax = "boot_CI_97.5",
                               color = point_color_col), size = 0.2) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(size = 1, fill = NA),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x.bottom = element_text(size = 10, color = "black"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 14, color = "black"),
      strip.text.x = element_text(angle = 0, color = "black"),
      strip.text.y = element_text(angle = 0, color = "black"),
    ) +
    scale_y_continuous(
      breaks = function(x) {
        auto_breaks <- pretty(x)
        extra_breaks <- c(-0.58, 0.58)
        sort(unique(c(auto_breaks, extra_breaks)))
      }) +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = point_color_palette) +
    labs(y = expression("Log2(FD)"), x = "", color = point_color_title)

  return(p)
}

grouped_permutation_test_plot <- function(
    sc_utils_obj,
    celltype_col,
    cell_group_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    panel_bg_cols = NULL,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5)
) {
  meta_data <- as.data.frame(sc_utils_obj@meta_data)
  permutation_data <- as.data.frame(sc_utils_obj@results$permutation)

  plot_data <- permutation_data %>%
    left_join(meta_data, by = c("clusters" = celltype_col)) %>%
    mutate(significance =
           ifelse(FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
              paste("FDR <", FDR_threshold, "&", 
                    expression("|Log2(FD)| >="),
                    round(log2FD_threshold, 2)), "n.s."),
           significance =
             factor(significance,
                    levels = c(paste("FDR <", FDR_threshold, "&",
                               expression("|Log2(FD)| >="),
                               round(log2FD_threshold, 2)), "n.s.")))

  plot_data$clusters <- factor(plot_data$clusters, levels = levels(meta_data[[celltype_col]]))

  bg_data <- unique(plot_data[, c(cell_group_col, "clusters")])
  bg_data$ymin <- -Inf
  bg_data$ymax <- Inf
  bg_data$obs_log2FD <- 0

  if (is.null(panel_bg_cols)) {
    unique_groups <- length(unique(plot_data[[cell_group_col]]))
    panel_bg_cols <- brewer.pal(min(unique_groups, 9), "Pastel1")
    names(panel_bg_cols) <- unique(plot_data[[cell_group_col]])
  } else {
    panel_bg_cols <- scales::alpha(panel_bg_cols, 0.2)
  }

  p <- ggplot(plot_data, aes_string(x = "clusters", y = "obs_log2FD")) +
    geom_rect(data = bg_data, aes_string(fill = cell_group_col),
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1) +
    geom_hline(yintercept = -1) +
    scale_fill_manual(values = panel_bg_cols, guide = 'none') +
    geom_pointrange(aes_string(ymin = "boot_CI_2.5", ymax = "boot_CI_97.5",
                               color = "significance")) +
    scale_y_continuous(breaks = c(-1, -0.58, 0, 0.58, 1)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(size = 1, fill = NA),
      axis.text.y = element_text(size = 11.5, color = "black"),
      axis.text.x.bottom = element_text(size = 11.5, color = "black"),
      strip.text = element_text(size = 14.5, color = "black"),
      strip.text.y = element_text(angle = 0, color = "black")
    ) +
    scale_color_manual(values = c("salmon", "grey")) +
    labs(y = expression("Log2(FD)"), x = "", color = "Significance") +
    facet_grid(as.formula(paste("~", cell_group_col)),
               scales = "free_x", space = "free")

  return(p)
}

apply_mods <- function(plt, mods) {
  if (!is.null(mods)) {
    if (is.list(mods)) {
      plt <- Reduce(`+`, mods, init = plt)
    } else {
      plt <- plt + mods
    }
  }
  return(plt)
}

combined_permutation_test_plot <- function(
    seurat_obj,
    prop_test,
    celltype_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    point_color_col = NULL,
    point_color_title = NULL,
    point_color_palette = NULL,
    perm_point_color_by_celltype = F,
    perm_point_color_title = NULL,
    perm_point_color_palette = NULL,
    coord_flip = F,
    barplot_log_scale = F,
    mods_barplot = NULL,
    mods_permplot = NULL
) {
  comparative_plot <-
    proportion_barplot(seurat_obj, celltype_col, condition_col,
                       condition_title, condition_cols, point_color_col,
                       point_color_title, point_color_palette,
                       barplot_log_scale) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  permutation_plot <-
    permutation_test_plot(prop_test, celltype_col, condition_col,
                          condition_title, condition_cols,
                          FDR_threshold, log2FD_threshold,
                          perm_point_color_by_celltype,
                          perm_point_color_title,
                          perm_point_color_palette)
  if (coord_flip) {
    comparative_plot <- comparative_plot + coord_flip()
    comparative_plot <- comparative_plot +
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.box.just = "left",
            legend.margin = margin(0, 0, 0, 0),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    permutation_plot <- permutation_plot + coord_flip()
    permutation_plot <- permutation_plot +
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.box.just = "left",
            legend.margin = margin(0, 0, 0, 0),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size = 10))
    comparative_plot <- apply_mods(comparative_plot, mods_barplot)
    permutation_plot <- apply_mods(permutation_plot, mods_permplot)
    combined_plot <- plot_grid(permutation_plot, comparative_plot,
                         nrow = 1, align = "h", axis = "tb")
  } else {
    comparative_plot <- comparative_plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    comparative_plot <- apply_mods(comparative_plot, mods_barplot)
    permutation_plot <- apply_mods(permutation_plot, mods_permplot)
    combined_plot <- plot_grid(comparative_plot, permutation_plot,
                               ncol = 1, align = "v", axis = "lr")
  }
  return(combined_plot)
}

grouped_combined_permutation_test_plot <- function(
    seurat_obj,
    prop_test,
    celltype_col,
    cell_group_col,
    condition_col,
    condition_title = condition_col,
    condition_cols = c("darkgreen", "darkred"),
    panel_bg_cols = NULL,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    point_color_col = NULL,
    point_color_title = NULL,
    point_color_palette = NULL,
    barplot_log_scale = F,
    mods_barplot = NULL,
    mods_permplot = NULL
) {
  comparative_plot <-
    grouped_proportion_barplot(seurat_obj, celltype_col, cell_group_col,
                               condition_col, condition_title, condition_cols,
                               panel_bg_cols, point_color_col,
                               point_color_title, point_color_palette,
                               barplot_log_scale) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

  permutation_plot <-
    grouped_permutation_test_plot(prop_test, celltype_col, cell_group_col,
                                  condition_col, condition_title,
                                  condition_cols, panel_bg_cols, FDR_threshold,
                                  log2FD_threshold) +
    theme(strip.text.x = element_blank())

  comparative_plot <- apply_mods(comparative_plot, mods_barplot)
  permutation_plot <- apply_mods(permutation_plot, mods_permplot)
  combined_plot <- plot_grid(comparative_plot, permutation_plot, ncol = 1,
                             align = "v", axis = "lr")

  return(combined_plot)
}
