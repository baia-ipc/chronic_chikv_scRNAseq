#
# (c) Giorgio Gonnella, Sopheap Oeng
#     Institut Pasteur du Cambodge, 2025
#     License: CC-BY-SA
#
# Wrappers to clusterProfiler ORA and GSEA visualization functions
#

library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(forcats)
library(stringr)

sanitize_name <- function(name) {
    cleaned_name <- gsub("[^a-zA-Z0-9_]", "_", name)
    cleaned_name <- gsub("_+", "_", cleaned_name)
    cleaned_name <- gsub("^_|_$", "", cleaned_name)
    return(cleaned_name)
}

wrap_description <- function(desc, lw) {
    words <- str_split(desc, " ")[[1]]
    num_words <- length(words)
    if (num_words <= lw) {
      return(desc)
    }
    result <- character()
    for (i in seq(1, num_words, by = lw)) {
      segment <- paste(words[i:min(i + lw - 1, num_words)], collapse = " ")
    result <- c(result, segment)
    }
    return(paste(result, collapse = "\n"))
}

gsea_allplots <- function(gsea_results, sorted_genes, dirname, filepfx,
                          show = T, lineplots = T) {
    try(gsea_barplot(gsea_results,
                 paste0(dirname, "/", filepfx, ".GSEA.barplot.pdf"),
                 show = show))
    try(gsea_dotplot(gsea_results,
                 paste0(dirname, "/", filepfx, ".GSEA.dotplot.pdf"),
                 show = show))
    try(gsea_cnetplot(gsea_results, sorted_genes,
                  paste0(dirname, "/", filepfx, ".GSEA.cnetplot.pdf"),
                  show = show))
    try(gsea_heatmap(gsea_results, sorted_genes,
                 paste0(dirname, "/", filepfx, ".GSEA.heatmap.pdf"),
                 show = show))
    try(gsea_emapplot(gsea_results,
                  paste0(dirname, "/", filepfx, ".GSEA.emapplot.pdf"),
                  show = show))
    try(gsea_ridgeplot(gsea_results,
                   paste0(dirname, "/", filepfx, ".GSEA.ridgeplot.pdf"),
                   show = show))
    if (lineplots) {
      gsea_lineplots(gsea_results,
                     paste0(dirname, "/", filepfx, ".GSEA.lineplots/"),
                     filepfx, show = show)
    }
}

gsea_barplot <- function(gsea_results, filename = NULL, n = 20,
                         colors = c("blue", "red"), show = T) {
    gse_df <- as.data.frame(gsea_results)
    top_terms <- gse_df[1:n, ]
    top_terms <- na.omit(top_terms)
    plt <- ggplot(top_terms, aes(x = reorder(Description, NES),
                                 y = NES, fill = p.adjust)) +
           geom_bar(stat = "identity") +
           coord_flip() +
           scale_fill_gradient(low = colors[[1]], high = colors[[2]]) +
           labs(title = paste0("Top ", n," Enriched Terms"),
                x = "Term", y = "Normalized Enrichment Score (NES)") +
           theme_minimal() +
           theme(axis.text.y = element_text(size = 10),
                 axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12))
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 10, height = 8)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

gsea_dotplot <- function(gsea_results, filename = NULL, n = 10, show = T) {
    plt <- dotplot(gsea_results, showCategory = n, split = ".sign") +
          facet_grid(.~.sign)
    if (!is.null(filename)) {
      pdf(filename, width = 8, height = 12)
      print(plt)
      dev.off()
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

gsea_cnetplot <- function(gsea_results, sorted_genes, filename = NULL,
                    n = 20, show = T) {
  plt <- cnetplot(gsea_results, categorySize = "pvalue",
              foldChange = sorted_genes,
              color_category = "firebrick",
              color_gene = "steelblue",
              showCategory = n, rescale = TRUE)
  if (!is.null(filename)) {
    pdf(filename, width = 8, height = 8)
    print(plt)
    dev.off()
  }
  if (show) {
    print(plt)
  }
  return(plt)
}

gsea_heatmap <- function(gsea_results, sorted_genes, filename = NULL,
                    n = 20, show = T) {
  plt <- heatplot(gsea_results, foldChange = sorted_genes,
              showCategory = n)
  if (!is.null(filename)) {
    pdf(filename, width = 10, height = 8)
    print(plt)
    dev.off()
  }
  if (show) {
    print(plt)
  }
  return(plt)
}

gsea_emapplot <- function(gsea_results, filename = NULL, n = 30, show = T) {
  plt <- emapplot(pairwise_termsim(gsea_results), showCategory = n)
  if (!is.null(filename)) {
    pdf(filename, width = 10, height = 8)
    print(plt)
    dev.off()
  }
  if (show) {
    print(plt)
  }
  return(plt)
}

gsea_ridgeplot <- function(gsea_results, filename = NULL, n = 20, show = T) {
  plt <- ridgeplot(gsea_results, showCategory = n) +
        labs(x = "Enrichment Distribution")
  if (!is.null(filename)) {
    pdf(filename, width = 10, height = 8)
    print(plt)
    dev.off()
  }
  if (show) {
    print(plt)
  }
  return(plt)
}

gsea_lineplots <- function(gsea_results, dirname = NULL, filepfx, show = T) {
  for (i in 1:nrow(gsea_results)) {
    plt <- try(gseaplot(gsea_results, by = "all",
                    title = gsea_results$Description[i], geneSetID = i))
    if (class(plt != "try-error")) {
      if (!is.null(dirname)) {
        cleaned_name <- sanitize_name(gsea_results$Description[i])
        ggsave(plot = plt, filename = paste0(dirname, "/", filepfx,
                                             ".gsea_plot.", cleaned_name,".pdf"),
               width = 10, height = 8, create.dir = TRUE)
      }
      if (show) {
          print(plt)
      }
    }
  }
}

ora_allplots <- function(ora_results, l2fc, dirname, filepfx,
                         show = T) {
    try(ora_barplot(ora_results,
                    paste0(dirname, "/", filepfx, ".ORA.barplot.pdf"),
                    show = show))
    try(ora_countsplot(ora_results,
                       paste0(dirname, "/", filepfx, ".ORA.countsplot.pdf"),
                       show = show))
    try(ora_qsbarplot(ora_results,
                      paste0(dirname, "/", filepfx, ".ORA.qsbarplot.pdf"),
                      show = show))
    try(ora_dotplot(ora_results,
                    paste0(dirname, "/", filepfx, ".ORA.dotplot.pdf"),
                    show = show))
    try(ora_cnetplot(ora_results,
                     paste0(dirname, "/", filepfx, ".ORA.cnetplot.pdf"),
                     show = show))
    try(ora_heatmap(ora_results, l2fc,
                    paste0(dirname, "/", filepfx, ".ORA.heatmap.pdf"),
                    show = show))
    try(ora_upsetplot(ora_results,
                      paste0(dirname, "/", filepfx, ".ORA.upsetplot.pdf"),
                      show = show))
    if (nrow(ora_results) >= 5) {
        try(ora_treeplot(ora_results,
                         paste0(dirname, "/", filepfx, ".ORA.treeplot.pdf"),
                         show = show))
    }
}

ora_barplot <- function(ora_results, filename = NULL, n = 20, show = T) {
    ora_results <- head(ora_results, n)
    plt <- ggplot(ora_results, aes(fct_reorder(Description, Count), Count)) +
            geom_bar(stat = "identity", fill = "skyblue") +
            coord_flip() +
            theme_bw() +
            theme(legend.position = "none") +
            labs(x = "Term", y = "Count",
                 title = paste0("Top ", n))
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 8, height = 6)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_countsplot <- function(ora_results, filename = NULL, n = 20, show = T) {
    ora_results <- as.data.frame(ora_results)
    edata <- ora_results[order(ora_results$p.adjust),]
    edata <- head(edata, n)
    plt <- ggplot(edata, aes(reorder(Description, Count), Count,
                                fill = p.adjust)) +
              geom_bar(stat = 'identity') +
              coord_flip() +
              scale_fill_gradient(low = "blue", high = "red") +
              theme_bw() +
              labs(x = "Term", y = "Counts",
                   title = paste0("Top ", n, " Enrichment"))
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 12, height = 10)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_qsbarplot <- function(ora_results, filename = NULL, n = 20, show = T) {
    ora_results <- as.data.frame(ora_results)
    edata <- ora_results[order(ora_results$p.adjust),]
    edata <- head(edata, n)
    plt <- mutate(edata, qscore = -log(p.adjust, base = 10)) %>%
               ggplot(aes(reorder(Description, qscore), qscore,
                          fill = p.adjust)) +
               geom_bar(stat = 'identity') +
               coord_flip() +
               scale_fill_gradient(low = "blue", high = "red") +
               theme_bw() +
               labs(x = "Terms", y = "qscore",
                    title = paste0("Top ", n, " enriched"))
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 12, height = 10)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_dotplot <- function(ora_results, filename = NULL, n = 20, show = T) {
    plt <- dotplot(ora_results, showCategory = n,
                   title = paste0("Top ", n, " enriched"))
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 10, height = 18)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_cnetplot <- function(ora_results, filename = NULL, n = 10, show = T) {
    plt <- cnetplot(ora_results, showCategory = n,
                       color_edge = "category")
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 15, height = 12)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_heatmap <- function(ora_results, l2fc, filename = NULL,
                    n = 15, show = T) {
    plt <- heatplot(ora_results, foldChange = l2fc,
                    showCategory = n)
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 10, height = 10)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_upsetplot <- function(ora_results, filename = NULL, n = 10, show = T) {
    plt <- upsetplot(ora_results, n = n)
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 15, height = 10)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

ora_treeplot <- function(ora_results, filename = NULL, show = T) {
    pwsim <- pairwise_termsim(ora_results)
    plt <- treeplot(pwsim)
    if (!is.null(filename)) {
        ggsave(filename, plot = plt, width = 14, height = 8)
    }
    if (show) {
        print(plt)
    }
    return(plt)
}

kegg_draw_pathways <- function(kegg_results, l2fc, dirname, verbose = T) {
    if (is.null(l2fc) || length(l2fc) == 0) {
        if (verbose) {
            print("Skipping KEGG pathway drawing: no gene fold-change data")
        }
        return(invisible(NULL))
    }
    sign_kegg <- kegg_results[kegg_results$p.adjust <= 0.05, ]
    pathway_ids <- rownames(sign_kegg)
    current_dir <- getwd()
    on.exit(setwd(current_dir), add = TRUE)
    for (id in pathway_ids) {
        row <- sign_kegg[sign_kegg$ID == id, ]
        desc <- sanitize_name(row$Description)
        if (verbose) {
            print(paste("Drawing pathway:", id))
        }
        setwd(dirname)
        try(pathview(gene.data = l2fc, pathway.id = id,
                     species = "hsa", kegg.native = FALSE, gene.idtype = "SYMBOL",
                     multi.state = F, out.suffix = paste0("KEGG.", desc),
                     kegg.dir = file.path(dirname)),
            silent = TRUE)
    }
    unlink(file.path(dirname, "*.xml"), recursive = TRUE)
}
