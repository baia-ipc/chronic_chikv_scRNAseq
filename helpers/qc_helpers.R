library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(reshape2)

#
# Note: the code of these functions is heavily based on code written by S.Mella
#

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

plot_higest_expressed <- function(so, nfeat2plot = 50,
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
