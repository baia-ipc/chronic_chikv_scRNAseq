#
# (c) Giorgio Gonnella, Sopheap Oeng
#     Institut Pasteur du Cambodge, 2025
#     License: CC-BY-SA
#
# Wrappers to clusterProfiler ORA and GSEA computation functions which:
#
# (1) use human as organism and gene symbols as identifiers
# (2) for ORA: select the significant genes as p-val adj < 0.05
# (3) for GSEA: sort the genes by log2FC for GSEA
# (4) map the gene symbols to Entrez IDs whenever necessary
#     reporting the number / percentage of symbols that cannot be translated
#

library(clusterProfiler)
library(dplyr)
library(pathview)
library(Organism.dplyr)
library(fgsea)
library(DOSE)
library(ReactomePA)
library(msigdbr)
require(org.Hs.eg.db)

### Computation parameters

MIN_GS_SIZE <- 3
MAX_GS_SIZE <- 800
P_VAL_CUTOFF <- 0.05
Q_VAL_CUTOFF <- 0.20
P_ADJ_METHOD <- "BH"
GSEA_EPS <- 1e-10
GSEA_EXPONENT <- 1

map_symbols_to_entrez_IDs <- function(symbols, verbose=T) {
    if (verbose) {
      print(paste0("Total number of genes: ", length(symbols)))
    }
    results <- c()
    if (length(symbols) > 0) {
      # try() necessary since mapids throws an exception if all are unmapped
      # possible alternative is bitr() of clusterProfiler, see below
      results <- try(mapIds(org.Hs.eg.db, keys = symbols,
                        column = "ENTREZID", keytype = "SYMBOL",
                        multiVals = "first"))
      if (class(results) == "try-error") {
        results = c()
      }
      if (verbose) {
        n_missing = sum(is.na(results))
        p_missing = round((n_missing / length(symbols)) * 100, 1)
        print(paste0("- of these, number of genes non considered, ",
                     "since missing Entrez ID: ", n_missing, " (",
                     p_missing, "%)"))
      }
    }
    return(results)
}

map_symbols_to_entrez_IDs_bitr <- function(symbols, verbose=T) {
    if (verbose) {
      print(paste0("Total number of genes: ", length(symbols)))
    }
    results <- c()
    if (length(symbols) > 0) {
      results <- bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db, drop=F)
      # Note: we need drop=F because we can use the retval to set names()
      if (verbose) {
        n_missing = sum(is.na(results))
        p_missing = round((n_missing / length(symbols)) * 100, 1)
        print(paste0("- of these, number of genes non considered, ",
                     "since missing Entrez ID: ", n_missing, " (",
                     p_missing, "%)"))
      }
    }
    return(results)
}

get_sorted_gene_symbols <- function(genes_info) {
    genes <- genes_info$avg_log2FC
    names(genes) <- rownames(genes_info)
    genes <- sort(genes, decreasing = TRUE)
    return(genes)
}

get_sorted_gene_entrez_IDs <- function(genes_info, verbose=T) {
    genes <- get_sorted_gene_symbols(genes_info)
    names(genes) <- map_symbols_to_entrez_IDs(names(genes), verbose)
    genes <- genes[!is.na(names(genes))]
    return(genes)
}

get_signif_gene_symbols <- function(genes_info) {
    genes_info <- subset(genes_info, p_val_adj <= 0.05)
    return(rownames(genes_info))
}

get_signif_gene_entrez_IDs <- function(genes_info, verbose=T) {
    genes_info <- subset(genes_info, p_val_adj <= 0.05)
    genes <- map_symbols_to_entrez_IDs(rownames(genes_info), verbose)
    genes <- genes[!is.na(genes)]
    return(genes)
}

get_MSigDb_term2gene <- function(category) {
    m_df <- msigdbr(species = "Homo sapiens", category = category)
    m_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)
    return(m_t2g)
}

get_signif_gene_info_symbols <- function(genes_info) {
    results <- subset(genes_info, p_val_adj <= 0.05)
    return(results)
}

get_signif_gene_info_entrez_IDs <- function(genes_info, verbose=T) {
    genes_info <- subset(genes_info, p_val_adj <= 0.05)
    entrez_IDs <- map_symbols_to_entrez_IDs(rownames(genes_info))
    results <- genes_info[!is.na(entrez_IDs),]
    rownames(results) <- entrez_IDs[!is.na(entrez_IDs)]
    return(results)
}

### Gene set enrichment analysis

# Arguments:
#   Data frame "genes_info", with a column avg_log2FC
#   where each row is a human gene (name: gene symbol)
#
# Return value:
#   Data frame containing the results of the GSEA analysis

compute_gsea_GO <- function(genes_info, verbose = T) {
    # verbose not used, but kept for consistency with the other functions
    genes <- get_sorted_gene_symbols(genes_info)
    results <- gseGO(geneList = genes,
                     ont = "ALL",
                     seed = TRUE,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     pvalueCutoff = 0.05,
                     minGSSize = 3,
                     maxGSSize = 800)
    return(results)
}

compute_gsea_KEGG <- function(genes_info, verbose = T) {
    genes <- get_sorted_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results <- gseKEGG(geneList = genes,
                       organism = "hsa",
                       seed = TRUE,
                       keyType = "ncbi-geneid",
                       pvalueCutoff = 0.05,
                       minGSSize = 3,
                       maxGSSize = 800)
    }
    return(results)
}

compute_gsea_Reactome <- function(genes_info, verbose = T) {
    genes <- get_sorted_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results <- gsePathway(geneList = genes,
                          organism = "human",
                          seed = TRUE,
                          exponent = 1,
                          minGSSize = 3,
                          maxGSSize = 800,
                          eps = 1e-10,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          verbose = verbose)
    }
    return(results)
}

compute_gsea_MSigDb <- function(genes_info, verbose = T, collection = NULL) {
    genes <- get_sorted_gene_symbols(genes_info)
    results <- GSEA(geneList = genes,
                    exponent = 1,
                    minGSSize = 3,
                    maxGSSize = 800,
                    eps = 1e-10,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = get_MSigDb_term2gene(collection),
                    verbose = verbose,
                    seed = TRUE,
                    by = "fgsea")
    return(results)
}
### Over-representation analysis

# Arguments:
#   Data frame "genes_info", with columns avg_log2FC and p_val_adj
#   where each row is a human gene (name: gene symbol)
#
# Return value:
#   Data frame containing the results of the GSEA analysis

compute_ora_GO <- function(genes_info, verbose = T) {
    genes <- get_signif_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results <- enrichGO(gene = genes,
                       universe = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       minGSSize = MIN_GS_SIZE,
                       maxGSSize = 800,
                       readable = TRUE)
    }
    return(results)
}

compute_ora_GGO_subset <- function(genes, ont, level) {
    ggo <- groupGO(genes, org.Hs.eg.db, ont = ont, level = level,
                   readable = TRUE)
    results <- ggo@result[order(-ggo@result$Count),]
    results$set <- ont
    results$level <- level
    return(results)
}

compute_ora_GGO <- function(genes_info, verbose = T) {
    genes <- get_signif_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results_BP2 <- compute_ora_GGO_subset(genes, "BP", 2)
      results_CC2 <- compute_ora_GGO_subset(genes, "CC", 2)
      results_MF2 <- compute_ora_GGO_subset(genes, "MF", 2)
      results_BP3 <- compute_ora_GGO_subset(genes, "BP", 3)
      results_CC3 <- compute_ora_GGO_subset(genes, "CC", 3)
      results_MF3 <- compute_ora_GGO_subset(genes, "MF", 3)
      results <- rbind(results_BP2, results_CC2, results_MF2,
                       results_BP3, results_CC3, results_MF3)
    }
    return(results)
}

compute_ora_KEGG <- function(genes_info, verbose = T) {
    genes_info <- genes_info[!is.na(genes_info$p_val_adj) &
                             !is.na(genes_info$avg_log2FC), ]
    genes <- get_signif_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results <- enrichKEGG(gene = genes,
                       universe = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       minGSSize = 3,
                       maxGSSize = 800)
    }
    return(results)
}

compute_ora_Reactome <- function(genes_info, verbose = T) {
    genes <- get_signif_gene_entrez_IDs(genes_info, verbose)
    results <- NULL
    if (length(genes) > 0) {
      results <- enrichPathway(gene = genes,
                       organism = "human",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       universe = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                       minGSSize = 3,
                       maxGSSize = 800,
                       readable = TRUE)
    }
    return(results)
}

compute_ora_MSigDb <- function(genes_info, verbose = T, collection = NULL) {
    # verbose not used, but kept for consistency with the other functions
    genes <- get_signif_gene_symbols(genes_info)
    results <- NULL
    if (length(genes) > 0) {
      results <- enricher(genes,
                        TERM2GENE=get_MSigDb_term2gene(collection),
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        minGSSize = 3,
                        maxGSSize = 800)
    }
    return(results)
}
