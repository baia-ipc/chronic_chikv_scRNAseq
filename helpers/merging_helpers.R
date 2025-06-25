so_version <- function(so) {
  so_v = substr(so@version, 1, 1)
  if (!so_v %in% c("4", "5"))
    stop("Unsupported Seurat version: ", so@version)
  return(so_v)
}

prepare_so_merging <- function(so) {
  # remove all assays except RNA
  DefaultAssay(so) <- "RNA"
  for (assay in names(so@assays)) {
    if (assay != "RNA") {
      so@assays[[assay]] <- NULL
    }
  }
  if (substr(so@version, 1, 1) == "5") {
    so@assays$RNA@layers$data <- NULL
    so@assays$RNA@layers$scale.data <- NULL
  }
  # remove SCT results from meta.data
  so@meta.data <-
    so@meta.data[-grep(pattern = "SCT", colnames(so@meta.data))]
  # remove all embeddings
  for (embedding in names(so@reductions)) {
    so@reductions[[embedding]] <- NULL
  }
  # remove all graphs
  for (graph in names(so@graphs)) {
    so@graphs[[graph]] <- NULL
  }
  return(so)
}

prepare_so_merging_Seurat4 <- function(so, verbose=F){
  if (verbose) cat("Changing default assay to 'RNA'\n")
  DefaultAssay(so) <- "RNA"
  if (verbose) cat("Removing SCT assay\n")
  so@assays$SCT <- NULL
  so@meta.data <-
    so@meta.data[-grep(pattern = "SCT", colnames(so@meta.data))]
  if (verbose) cat("Removing graphs\n")
  for(ng in names(so@graphs)){
    so[[ng]] <- NULL
  }
  for(r in Seurat::Reductions(so)){
    so[[r]] <- NULL
  }
  if (verbose) cat("\t",unique(so$orig.ident, " done.\n"))
  return(so)
}
