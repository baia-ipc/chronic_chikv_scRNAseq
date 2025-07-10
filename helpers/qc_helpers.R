library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(reshape2)

#'
#' Quality control functions.
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023-2025

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

test_cell_cycle_effect <- function(so, so_label, s_genes, g2m_genes,
                                   verbose = T, showplots = F) {

  if (verbose) cat("\n### Processing ", so_label, "\n")

  if (verbose) cat("\n##### Number of S genes in ", so_label, " : \n")
  s_genes_table <- as.data.frame(table(s_genes %in% rownames(so)))
  if (verbose)
    cat("\n", paste(apply(s_genes_table,1, paste, collapse = " : "),
                    collapse  = "\n\n"))

  if (verbose) cat("\n\n##### Number of G2M genes in ", so_label, " : \n")
  g2m_genes_table <- as.data.frame(table(g2m_genes %in% rownames(so)))
  if (verbose)
    cat("\n", paste(apply(g2m_genes_table,1, paste, collapse = " : "),
                    collapse  = "\n\n"))

  if (verbose) cat("\n\n#### Normalizing Seurat Object by nUMI count\n")
  so <- NormalizeData(so, assay = "RNA")

  if (verbose) cat("\n#### Score cells for cell cycle\n")
  so <- CellCycleScoring(object = so,
                         g2m.features = g2m_genes[g2m_genes %in% rownames(so)],
                         s.features = s_genes[s_genes %in% rownames(so)])

  if (verbose) cat("\n##### Frequency of the different cell cycle phases\n")
  ccp_df <- as.data.frame.array(table(so$Phase))
  ccp_df$phase <- rownames(ccp_df)

  if (verbose) cat("\n", paste(apply(ccp_df[,2:1],1, paste, collapse = " : "),
                               collapse  = "\n\n"))

  if (verbose) cat("\n\n#### Preparing data for PCA\n")
  if (verbose) cat("\nIdentify most variable genes\n")
  so <- FindVariableFeatures(so, selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE, assay = "RNA")

  if (verbose) cat("\n#### Scaling the counts\n")
  so <- ScaleData(so, assay = "RNA", verbose = FALSE)

  if (verbose) cat("\n#### Perform PCA\n")
  so <- RunPCA(so, assay = "RNA", reduction.name = "pcaRNA",
               reduction.key = "pcaRNA_", verbose = FALSE)

  if (verbose) cat("\n#### Plot the PCA colored by cell cycle phase\n")
  ccg_colors <-
    RColorBrewer::brewer.pal(length(levels(so$Phase)), name = "Set1")

  if (showplots) {
    plot(
      DimPlot(so, reduction = "pcaRNA", group.by= "Phase", cols = ccg_colors)
    )
    plot(
      DimPlot(so, reduction = "pcaRNA", group.by= "Phase", split.by = "Phase",
              cols = ccg_colors)
    )
  }
  return(so)
}

compute_percent_mito <- function(so, mitoGenes){
  percent.mito <-
    Matrix::colSums(so@assays$RNA$counts[
                      rownames(so@assays$RNA$counts) %in% mitoGenes,]) /
                    Matrix::colSums(so@assays$RNA$counts)
  so@meta.data$percent.mito <- percent.mito
  so@meta.data$percent.mito.cl <- cut(so@meta.data$percent.mito,
                                      breaks = c(-1, .1, .15, .2, 1))
  return(so)
}

plot_basic_qc <- function(so, title,
                          colors = c("green1", "blue1", "orange1","red1"),
                          alpha = .6) {
  p.pMito <- ggplot(so@meta.data,
                    aes(x = levels(factor(orig.ident)), y = percent.mito)) +
    geom_jitter(aes(color = percent.mito.cl), alpha = alpha) +
    geom_violin(trim = T, width = .5, alpha = .5) +
    scale_color_manual(values = colors) +
    xlab("") +
    guides(color = F) +
    theme_minimal()

  p.nGene <- ggplot(so@meta.data,
                    aes(x = levels(factor(orig.ident)), y = nFeature_RNA)) +
    geom_jitter(aes(color = percent.mito.cl), alpha = alpha) +
    geom_violin(trim = T, width = .5, alpha = .5) +
    scale_color_manual(values = colors) +
    xlab("") +
    guides(color = F) +
    ggtitle(title) +
    theme_minimal()

  p.nUMI <- ggplot(so@meta.data,
                   aes(x = levels(factor(orig.ident)), y = nCount_RNA)) +
    geom_jitter(aes(color = percent.mito.cl), alpha = alpha) +
    geom_violin(trim = T, width = .5, alpha = .5) +
    scale_color_manual(values = colors) +
    xlab("") +
    guides(color = F) +
    theme_minimal()

  grid.arrange(p.pMito, p.nGene, p.nUMI, ncol = 3, nrow= 1)
}

factor_resolution <- function(so){
  tmp_df <- so@meta.data
  ind <- grep("res.", colnames(tmp_df))
  for(i in ind){
    tmp_df[,i] <- factor(tmp_df[,i])
    tmp_df[,i] <- factor(tmp_df[,i],
                         levels = 0 : (length(levels(tmp_df[,i])) - 1) )
  }
  so@meta.data <- tmp_df
  return(so)
}

plot_highest_expressed <- function(so, nfeat2plot = 50,
                                  boxplot_alpha = 0.75,
                                  ylab_txt = "counts",
                                  colors = colorRampPalette(
                                        rev(pal_futurama("planetexpress"
                                                         )(12)))(nfeat2plot)) {
  mat_exp <- as.matrix(so@assays$RNA$counts)

  # calculate sum of expression (sum of UMI)
  ave_exprs <- rowSums2(mat_exp)

  # ordering vector according to the level of detected genes
  ord_feat <- order(ave_exprs, decreasing = TRUE)
  selfeatures <- head(ord_feat, nfeat2plot)

  # Subtracting the expression matrix
  sub_mat <- as.data.frame(mat_exp[selfeatures, ])
  sub_mat$feature <- factor(rownames(sub_mat), levels = rev(rownames(sub_mat)))

  # converting dataframe into long format
  sub_mat_l <- melt(sub_mat, id.vars = "feature")

  # make the plot
  plt <- ggplot(sub_mat_l, aes(x = value, y = feature)) +
    geom_boxplot(aes(fill = feature), alpha = boxplot_alpha) +
    scale_fill_manual(values = colors) +
    theme_light() +
    guides(fill = FALSE) +
    ylab(ylab_txt)
  plt
}

#' Add QC information from scater library to Seurat object
#'
#' @param so Seurat object
#'
#' @return Seurat object with additional QC metrics
scater_qc <- function(so) {
  sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = so@assays$RNA@counts))
  is_mito <- grep("^MT-", rownames(sce))
  qcstats <- scater::perCellQCMetrics(sce, subsets=list(Mito=is_mito))
  filt_ad <- data.frame(
    qc.lib.low = scater::isOutlier(qcstats$sum, log = TRUE, type = "lower"),
    qc.lib.high = scater::isOutlier(qcstats$sum,
                            log = FALSE, type = "higher", nmads = 4),
    qc.nexpr.low = scater::isOutlier(qcstats$detected, log = TRUE,
                                     type = "lower"),
    qc.nexpr.high = scater::isOutlier(qcstats$detected,
                              log = FALSE, type = "higher", nmads = 4),
    qc.mito.high = scater::isOutlier(qcstats$subsets_Mito_percent,
                             log = FALSE, type = "higher")
  )
  qcstats <- cbind(qcstats, filt_ad)
  names(qcstats) <- paste0("scater_", names(qcstats))
  so@meta.data <- cbind(so@meta.data, qcstats)
  withr::local_seed(1234)
  stats <- cbind(log10(so$nCount_RNA), log10(so$nFeature_RNA), so$mitoRatio)
  outlying <- robustbase::adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- scater::isOutlier(outlying, type = "higher")
  so@meta.data <- cbind(so@meta.data, multi.outlier)
  so
}

#' QC metric distribution histogram
#'
#' An histogram plot of the distribution
#' for one of the QC metrics in a Seurat object
#'
#' @param so Seurat object
#' @param measure String. The name of the metric to plot
#'
#' @return A ggplot object
#'
qc_distri_plot <- function(so, measure) {
  ggplot2::ggplot(so@meta.data, ggplot2::aes_string(x=measure)) +
    ggplot2::geom_histogram(bins=100) +
  ggplot2::ggtitle(paste(measure, "Distribution")) +
  ggplot2::xlab(measure) + ggplot2::ylab("Frequency")
}


#' Plot number of cells in different samples
#'
#' @param so Seurat object
#'
histogram_n_cells <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=sample, fill=sample)) +
    	ggplot2::geom_bar() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
                     ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of cells")
}

#' Plot number of UMIs per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_umis <- function(so, plot_nrows=0) {
  so@meta.data %>%
      ggplot2::ggplot(ggplot2::aes(color=sample, x=nCount_RNA, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::scale_x_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::ylab("cell density") +
    	ggplot2::geom_vline(xintercept = 100) +
    	ggplot2::geom_vline(xintercept = 500) +
    	ggplot2::geom_vline(xintercept = 1000) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
    	ggplot2::ggtitle("Number of UMIs per cell")
}

#' Boxplot of the number of UMIs per cell
#'
#' @param so Seurat object
#'
boxplot_n_umis <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(
          ggplot2::aes(x=sample, y=log10(nCount_RNA), fill=sample)) +
    	ggplot2::geom_boxplot() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
          ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
          ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of UMIs per cell")
}

#' Plot number of genes per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_genes <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(color=sample, x=nFeature_RNA, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::theme_classic() +
    	ggplot2::scale_x_log10() +
    	ggplot2::geom_vline(xintercept = 300) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
    	ggplot2::ggtitle("Number of genes per cell")
}


#' Boxplot of the number of genes per cell
#'
#' @param so Seurat object
#'
boxplot_n_genes <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=sample,
                                   y=log10(nFeature_RNA), fill=sample)) +
    	ggplot2::geom_boxplot() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
        ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
        ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of genes per cell")
}

dotplot_n_umis_genes_mito <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=nCount_RNA,
                                   y=nFeature_RNA, color=mitoRatio)) +
    	ggplot2::geom_point() +
  	  ggplot2::scale_colour_gradient(low = "gray90", high = "black") +
    	ggplot2::stat_smooth(method=lm) +
    	ggplot2::scale_x_log10() +
    	ggplot2::scale_y_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 500) +
    	ggplot2::geom_hline(yintercept = 250) +
    	ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("N. UMIs vs N. genes and mit.genes ratio")
}

density_plot_mito_ratio <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(color=sample, x=mitoRatio, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::scale_x_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 0.2) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("Mitochondrial gene ratio per cell")
}

density_plot_complexity <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(
        ggplot2::aes(x=log10GenesPerUMI, color=sample, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 0.8) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("Transcriptional complexity")
}

#' Create and print several QC metrics plots
#'
#' Creates several QC metrics plots for a Seurat object
#' and prints them to the screen
#'
#' @param so Seurat object
#'
show_qc_plots <- function(so, plot_nrows=0) {
  print(histogram_n_cells(so))
  print(density_plot_n_umis(so, plot_nrows))
  print(boxplot_n_umis(so))
  print(density_plot_n_genes(so, plot_nrows))
  print(boxplot_n_genes(so))
  print(density_plot_mito_ratio(so, plot_nrows))
  print(dotplot_n_umis_genes_mito(so, plot_nrows))
  print(density_plot_complexity(so, plot_nrows))
}

