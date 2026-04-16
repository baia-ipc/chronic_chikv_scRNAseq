library(ggplot2)

# Return the two condition labels from a comparison label such as "A.NC_vs_C"
# or "NC_vs_C".
de_conditions_from_pair_label <- function(pair_label) {
  contrast <- sub("^[^.]+\\.", "", pair_label)
  conditions <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  if (length(conditions) != 2L) {
    stop("Cannot parse compared conditions from pair label: ", pair_label)
  }
  conditions
}

# Count the number of rows with positive and negative log2 fold-change per
# cell type. Pass a p-value-filtered list to count significant DEGs only.
de_direction_count_data <- function(de_results, celltypes = NULL,
                                    log2fc_col = "avg_log2FC",
                                    positive_label = "Up",
                                    negative_label = "Down",
                                    drop_empty = TRUE,
                                    sort_by_total = TRUE) {
  if (!is.null(de_results$pval_filt)) {
    de_results <- de_results$pval_filt
  }
  if (is.null(celltypes)) {
    celltypes <- names(de_results)
  }
  if (length(celltypes) == 0L) {
    return(data.frame(celltype = character(0), direction = character(0),
                      n = integer(0), total = integer(0)))
  }

  rows <- lapply(celltypes, function(ct) {
    ct_data <- de_results[[ct]]
    if (is.null(ct_data) || nrow(ct_data) == 0L || !log2fc_col %in% colnames(ct_data)) {
      pos_n <- 0L
      neg_n <- 0L
    } else {
      log2fc <- ct_data[[log2fc_col]]
      pos_n <- sum(!is.na(log2fc) & log2fc > 0)
      neg_n <- sum(!is.na(log2fc) & log2fc < 0)
    }
    data.frame(
      celltype = ct,
      direction = c(positive_label, negative_label),
      n = c(pos_n, neg_n),
      stringsAsFactors = FALSE
    )
  })

  counts <- do.call(rbind, rows)
  totals <- aggregate(n ~ celltype, counts, sum)
  names(totals)[[2]] <- "total"
  counts <- merge(counts, totals, by = "celltype", all.x = TRUE, sort = FALSE)

  if (drop_empty) {
    counts <- counts[counts$total > 0, , drop = FALSE]
  }
  if (nrow(counts) == 0L) {
    counts$celltype <- factor(counts$celltype)
    counts$direction <- factor(counts$direction,
                               levels = c(positive_label, negative_label))
    return(counts)
  }

  if (sort_by_total) {
    celltype_levels <- totals$celltype[order(totals$total, decreasing = TRUE)]
  } else {
    celltype_levels <- celltypes
  }
  celltype_levels <- celltype_levels[celltype_levels %in% unique(counts$celltype)]
  counts$celltype <- factor(counts$celltype, levels = celltype_levels)
  counts$direction <- factor(counts$direction,
                             levels = c(positive_label, negative_label))
  counts
}

de_direction_count_plot <- function(de_results, positive_condition,
                                    negative_condition, celltypes = NULL,
                                    title = NULL,
                                    log2fc_col = "avg_log2FC",
                                    y_label = "Number of DEGs",
                                    colors = c("#6BAED6", "#C44E52"),
                                    drop_empty = TRUE,
                                    sort_by_total = TRUE) {
  positive_label <- paste0("Up in ", positive_condition)
  negative_label <- paste0("Up in ", negative_condition)
  counts <- de_direction_count_data(
    de_results,
    celltypes = celltypes,
    log2fc_col = log2fc_col,
    positive_label = positive_label,
    negative_label = negative_label,
    drop_empty = drop_empty,
    sort_by_total = sort_by_total
  )

  if (nrow(counts) == 0L) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No significant DEGs") +
        theme_void() +
        labs(title = title)
    )
  }

  ggplot(counts, aes(x = celltype, y = n, fill = direction)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(
      values = stats::setNames(colors, c(positive_label, negative_label)),
      breaks = c(positive_label, negative_label)
    ) +
    labs(title = title, x = NULL, y = y_label, fill = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    )
}
