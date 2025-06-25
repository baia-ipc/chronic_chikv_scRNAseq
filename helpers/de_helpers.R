library(ggplot2)
library(gridExtra)
library(pheatmap)

run_pb_deseq2 <- function(pseudo, selgenes, ident.1, ident.2,
                          freqs = NULL) {
  print(paste0("Computing DESeq2 for ", ident.1, " vs ", ident.2, "..."))
  if (!ident.1 %in% Idents(pseudo)) {
    print(paste0("  Skipped as ", ident.1 ," not found in pseudo"))
    return()
  }
  if (!ident.2 %in% Idents(pseudo)) {
    print(paste0("  Skipped as ", ident.2 ," not found in pseudo"))
    return()
  }
  if (is.null(freqs)) {
    freqs <- table(Idents(pseudo))
  }
  print(paste0("  Frequency of ", ident.1, ": ", freqs[[ident.1]]))
  print(paste0("  Frequency of ", ident.2, ": ", freqs[[ident.2]]))
  if (freqs[[ident.1]] < 3) {
    print(paste0("  Skipped as ", ident.1 ," has less than 3 cells"))
    return()
  }
  if (freqs[[ident.2]] < 3) {
    print(paste0("  Skipped as ", ident.2 ," has less than 3 cells"))
    return()
  }

  markers <- FindMarkers(pseudo, features = selgenes, logfc.threshold = 0.1,
                         ident.1 = ident.1, ident.2 = ident.2,
                         test.use = "DESeq2")
  results <- list()
  results$all <- markers
  results$pval_filt <- markers %>% filter(p_val_adj < 0.05)
  print(paste0("  Number of genes with adj p-value < 0.05: ",
               nrow(results$pval_filt)))
  results
}

#
# Create a heatmap for a list of genes in multiple cell types,
# starting from the output of the DE analysis.
#
# Mandatory arguments:
#   markers: a named list of data frames; the names shall be the cell types
#            the data frames shall have at least a column called "gene" and
#            a column named "avg_log2FC"
#   genes_of_interest: a list of genes to be visualized
#   output_file: name of the file to which to output the plot
#
# Optional arguments:
#   cell_types: a vector of cell types to be visualized; the default is all
#               cell types contained in markers; if the vector is named,
#               the names are used in the plot X axis (otherwise the values
#               are used)
#   max_abs_value: maximum value of the scala (and -max_abs_value will be the
#                  minimum); default is NULL, i.e. the value will be set
#                  automatically based on the data
#   color1, color2, color3: min, middle and max points of the color gradient
#   title: main label of the plot
#   add_stars: if to add stars for the adj_p_val significance [default: T]
#
DE_selected_genes_heatmap <- function(markers,
                                      genes_of_interest,
                                      output_file,
                                      cell_types = NULL,
                                      max_abs_value = NULL,
                                      color1 = "red",
                                      color2 = rgb(240/255, 240/255, 240/255),
                                      color3 = "blue",
                                      title = "Heatmap of selected genes",
                                      add_stars = T) {
  
  heatmap_data <- data.frame(gene = genes_of_interest)
  
  if (is.null(cell_types)) {
    cell_types = names(markers)
  }
  
  if (is.null(names(cell_types))) {
    names(cell_types) <- cell_types
  }
  
  for (cell_type_name in names(cell_types)) {
    cell_type = cell_types[[cell_type_name]]
    cell_data <- markers[[cell_type]]
    cell_data$gene <- rownames(cell_data)
    filtered_data <- cell_data[cell_data$gene %in% genes_of_interest, ]
    filtered_data <- filtered_data[match(genes_of_interest, filtered_data$gene), ]

    # Add stars based on p_val_adj
    log2fc_with_stars <- sapply(1:nrow(filtered_data), function(i) {
      log2fc <- filtered_data$avg_log2FC[i]
      p_val_adj <- filtered_data$p_val_adj[i]
      if (!add_stars) {
        return(paste0(log2fc, "| "))
      }
      if (is.na(p_val_adj)) {
        return(paste0(log2fc, "| "))
      }
      if (p_val_adj < 0.001) {
        return(paste0(log2fc, "|***"))
      } else if (p_val_adj < 0.01) {
        return(paste0(log2fc, "|**"))
      } else if (p_val_adj < 0.05) {
        return(paste0(log2fc, "|*"))
      } else {
        return(paste0(log2fc, "| "))
      }
    })
    
    heatmap_data[[cell_type_name]] <- log2fc_with_stars
  }

  rownames(heatmap_data) <- heatmap_data$gene
  
  heatmap_data <- heatmap_data[, -1]
  
  numeric_heatmap_data <- apply(heatmap_data, 2, function(x) as.numeric(sapply(strsplit(x, "\\|"), `[`, 1)))
  rownames(numeric_heatmap_data) <- rownames(heatmap_data)
  
  annotation_matrix <- apply(heatmap_data, 2, function(x) sapply(strsplit(x, "\\|"), `[`, 2))
  
  if (is.null(max_abs_value)) {
    max_abs_value <- max(abs(numeric_heatmap_data), na.rm = TRUE)
  }
  breaks <- seq(-max_abs_value, max_abs_value, length.out = 100)
  
  pheatmap(
    numeric_heatmap_data,
    scale = "none",
    cluster_rows = F,
    cluster_cols = F,
    color = colorRampPalette(c(color1, color2, color3))(length(breaks) - 1),
    breaks = breaks,
    main = title,
    display_numbers = annotation_matrix,
    filename = output_file,
    number_color = "black",
    fontsize_number = 14
  )
}

#
# Compute the maximum values of the scaled expression and of the percentage of
# expressing cells for a list of genes in a Seurat object.
#
# Arguments:
#   so: Seurat object
#   genes: list of genes to be considered
#   group.by: the variable to group the data by
#
max_expression <- function(so, genes, group.by) {
  to_plot_f <- SeuratObject::FetchData(so, vars = c(group.by, genes))
  to_plot_f <- tidyr::pivot_longer(to_plot_f, cols = c(2:ncol(to_plot_f)))
  to_plot_f <- dplyr::group_by(to_plot_f, !!as.name(group.by), name) %>%
    dplyr::summarise(mean_scaled_expression = mean(value[value > 0]),
                     perc_expressing_cells =
                       (length(value[value > 0])/dplyr::n() * 100))
  list(
    max_expr = max(to_plot_f$mean_scaled_expression, na.rm = TRUE),
    max_perc = max(to_plot_f$perc_expressing_cells, na.rm = TRUE)
  )
}

GEX_dotplot <- function(so, condcol, genes, cond1 = NULL, cond2 = NULL,
                        max_expr = 1, max_perc = 100, threshold.to.plot = 5,
                        title = "", remove_empty = T) {
  if (is.null(cond1)) {
    cond1 = levels(so@meta.data[[condcol]])[[1]]
  }
  if (is.null(cond2)) {
    cond2 = levels(so@meta.data[[condcol]])[[2]]
  }

  cells_1 = which(so@meta.data[[condcol]] == cond1)
  if (length(cells_1) == 0) {
    stop(paste0("No cells found for condition: ", cond1))
  }
  cells_2 = which(so@meta.data[[condcol]] == cond2)
  if (length(cells_2) == 0) {
    stop(paste0("No cells found for condition: ", cond1))
  }
  GEX_1 = subset(so, cells = cells_1)
  GEX_2 = subset(so, cells = cells_2)
  group <- NULL
  name <- NULL
  value <- NULL
  mean_scaled_expression <- NULL
  perc_expressing_cells <- NULL

  to_plot_1 <- SeuratObject::FetchData(GEX_1, vars = c(condcol, genes))
  to_plot_1$condition <- cond1

  to_plot_2 <- SeuratObject::FetchData(GEX_2, vars = c(condcol, genes))
  to_plot_2$condition <- cond2

  to_plot_f <- rbind(to_plot_1, to_plot_2)
  to_plot_f <- tidyr::pivot_longer(to_plot_f, cols = c(2:(ncol(to_plot_f) - 1)))
  names(to_plot_f)[1] <- "group"
  names(to_plot_f)[2] <- "condition"

  to_plot_sum_f <- to_plot_f %>%
    dplyr::group_by(group, condition, name) %>%
    dplyr::summarise(mean_scaled_expression = mean(value[value > 0]),
                     perc_expressing_cells = (length(value[value > 0])/dplyr::n()*100))

  to_plot_sum_f$name <- ordered(as.factor(to_plot_sum_f$name), levels = genes)
  to_plot_sum_f$condition <- factor(to_plot_sum_f$condition, levels = c(cond1, cond2))

  to_plot_sum_f$perc_expressing_cells[to_plot_sum_f$perc_expressing_cells < threshold.to.plot] <- NA

  # remove data where both conditions have NA in perc_expressing_cells
  if (remove_empty) {
    to_plot_sum_f <- to_plot_sum_f[!is.na(to_plot_sum_f$perc_expressing_cells), ]
  }

  plot_out <- ggplot(to_plot_sum_f, aes(x = name, y = condition,
                                        col = mean_scaled_expression,
                                        size = perc_expressing_cells)) +
              geom_point(show.legend = T, position = position_jitter(width = 0, height = 0)) +
              cowplot::theme_cowplot() +
              theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                    legend.position = "bottom", legend.direction = "horizontal",
                    axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1),
                    plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
                    legend.box = "vertical") +
              labs(title = paste0("Expression by ", condcol), x = NULL,
                   y = NULL, color = "Mean scaled expression",
                   size = "% of expressing cells") +
              scale_color_viridis_c(option = "B", begin = 0.1, end = 0.9,
                                    limits = c(0, max_expr)) +
              scale_size_binned(range = c(1,9.5), limits = c(0, max_perc)) +
              scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
              ggtitle(title)

  plot_width <- 2 + length(genes)*0.4
  plot_height <- 2.4 + 0.15*max(nchar(genes))

  final_plot <- grid.arrange(plot_out, ncol=1,
                             widths = unit(plot_width, "in"),
                             heights = unit(plot_height, "in"))

  return(final_plot)
}

#
# Compare the expression of split.by column values [[1]] vs [[2]]
# in each of the groups in the group.by column.
#
cmpexpr_by_group <- function(so, genes_of_interest, split.by, group.by) {
  missing_genes <- setdiff(genes_of_interest, rownames(so))
  if (length(missing_genes) > 0) {
    stop(paste("The following genes are not found in the Seurat object:",
               paste(missing_genes, collapse = ", ")))
  }

  raw_counts <- GetAssayData(so, slot = "data")
  metadata <- so@meta.data

  if (!(split.by %in% colnames(metadata))) {
    stop(paste("Split column", split.by, "not found in metadata"))
  }

  if (!(group.by %in% colnames(metadata))) {
    stop(paste("Group column", group.by, "not found in metadata"))
  }

  results_list <- list()
  for (gene in genes_of_interest) {
    df <- data.frame(
      expression = raw_counts[gene, ],
      split = metadata[[split.by]],
      group = metadata[[group.by]]
    )
    mean_expression <- aggregate(expression ~ group + split,
                                 data = df, FUN = mean)
    mean_expression_wide <- reshape(mean_expression, idvar = "group",
                                    timevar = "split", direction = "wide")

    if (ncol(mean_expression_wide) > 2) {
      splits <- colnames(mean_expression_wide)[-1]
      avg_log2fc <- log2(mean_expression_wide[, splits[1]] /
                         mean_expression_wide[, splits[2]])
      mean_expression_wide$avg_log2fc <- avg_log2fc
    } else {
      mean_expression_wide$avg_log2fc <- NA
    }

    results_list[[gene]] <- mean_expression_wide
  }

  result_df <- do.call(rbind, lapply(names(results_list), function(gene) {
    df <- results_list[[gene]]
    colnames(df) <- gsub("expression.", "", colnames(df))
    cbind(Gene = gene, df)
  }))

  return(result_df)
}


#
# Barplot from the GSEA results
#
# Arguments:
#   gse: results of gseGO or gseKEGG
#
# Optional arguments:
#   top_n: number of top terms to plot [default: 20]
#   colors: vector of two colors for the gradient [default: blue, red]
#
gsea_barplot <- function(gse, top_n = 20, colors = c("blue", "red")) {
  gse_df <- as.data.frame(gse)
  top_terms <- gse_df[1:top_n, ]
  plt <- ggplot(top_terms, aes(x = reorder(Description, NES),
                             y = NES, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = colors[[1]], high = colors[[2]]) +
    labs(title = paste0("Top ", top_n," Enriched Terms"),
                 x = "Term", y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  return(plt)
}

#
# Select genes to use for analyses, such as DE analysis.
#
# Arguments:
#  so: Seurat object
#
# Optional arguments:
#   minrowsum: minimum row sum to keep a gene [default: 0]
#   keep_TCR: keep T cell receptor genes [default: F]
#   keep_IG: keep Immunoglobulin genes [default: F]
#   keep_MT: keep Mitocondrial genes [default: F]
#   keep_RP: keep Ribosomal protein genes [default: F]
#   keep_MALAT1: keep MALAT1 gene [default: F]
#
select_genes <- function(so, minrowsum = 0, keep_TCR = F,
                         keep_IG = F, keep_MT = F, keep_RP = F,
                         keep_MALAT1 = F, verbose = F) {
  selgenes <- rownames(so)
  if (verbose) {
    print(paste0("Initial number of genes: ", length(selgenes)))
  }
  if (minrowsum > 0) {
    if (verbose) {
      print(paste0("Minimum row sum parameter: ", minrowsum))
    }
    rs <- Matrix::rowSums(so@assays$RNA$counts)
    selgenes <- selgenes[rs > minrowsum]
    if (verbose) {
      print(paste0("Number of selected genes passing rowsum filter: ",
                   length(selgenes)))
    }
  }
  if (!keep_TCR) {
    TCRgenes <- grep(pattern = "^TR[AB][VC]", selgenes, value = TRUE)
    selgenes <- selgenes[!selgenes %in% TCRgenes]
    if (verbose) {
      print(paste0("Number of TCR genes to filter out: ", length(TCRgenes)))
      print(paste0("Number of selected genes after TCR filter: ",
                   length(selgenes)))
    }
  }
  if (!keep_IG) {
    IGgenes <- c(grep(pattern = "^IG[HLK][VDJ]", selgenes, value = TRUE),
                   c("IGHG1", "IGHD", "IGHE", "IGHA[12]", "IGHG[1234]",
                     "IGKC", "IGLC[1234567]", "AC233755.1"))
    selgenes <- selgenes[!selgenes %in% IGgenes]
    if (verbose) {
      print(paste0("Number of IG genes to filter out: ", length(IGgenes)))
      print(paste0("Number of selected genes after IG filter: ",
                   length(selgenes)))
    }
  }
  if (!keep_MT) {
    MTgenes <- grep(pattern = "^MT-", selgenes, value = TRUE)
    selgenes <- selgenes[!selgenes %in% MTgenes]
    if (verbose) {
      print(paste0("Number of Mitochondrial genes to filter out: ",
                   length(MTgenes)))
      print(paste0("Number of selected genes after MT filter: ",
                   length(selgenes)))
    }
  }
  if (!keep_RP) {
    RPgenes <- grep(pattern = "^RP[LS]", selgenes, value = TRUE)
    selgenes <- selgenes[!selgenes %in% RPgenes]
    if (verbose) {
      print(paste0("Number of Ribosomal protein genes to filter out: ",
                   length(RPgenes)))
      print(paste0("Number of selected genes after RP filter: ",
                   length(selgenes)))
    }
  }
  if (!keep_MALAT1) {
    selgenes <- selgenes[-which(selgenes=="MALAT1")]
    if (verbose) {
      print(paste0("Number of MALAT1 genes to filter out: 1"))
      print(paste0("Number of selected genes after MALAT1 filter: ",
                   length(selgenes)))
    }
  }
  if (verbose) {
    print(paste0("Total number of genes filtered out: ",
                 length(rownames(so)) - length(selgenes)))
    print(paste0("Final number of genes: ", length(selgenes)))
  }
  return(selgenes)
}
