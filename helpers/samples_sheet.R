#'
#' Handling of the samples sheet data
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#' Read samples sheet
#'
#' The sample sheet must be in TSV format. The first line must be
#' the header, containing column names and optionally starting with
#' a '#' (without any spaces following it).
#'
#' @param samples_sheet String. Path to the TSV file
#' @param column String. Column containing the sample IDs
#' @param verbose Boolean. Whether to print the sample sheet and the IDs.
#'
#' @return data.frame. The sample sheet.
#'
#' @examples
#'   samples_sheet <- read_samples_sheet("~/project/samples.tsv", "PatientID")
#'
read_samples_sheet <- function(filename, column="SampleID", verbose=TRUE) {
  header <- read.table(samples_sheet, sep="\t", header=FALSE,  nrows=1,
                       stringsAsFactors=FALSE, comment.char="")
  header <- gsub("^#", "", unlist(header))
  samples_metadata <- read.table(samples_sheet, sep="\t", col.names = header,
                                 header=FALSE, skip=1, comment.char="",
                                 colClasses = setNames("character", column))
  if (verbose) {
    print(paste("Sample sheet:", samples_sheet))
    print(samples_metadata)
  }

  samples_metadata
}

#' Get sample IDs from a sample sheet
#'
#' The sample sheet can be provided as a data frame, or as a path to a TSV file
#' (see the documentation of read_samples_sheet for details of the format).
#'
#' @param samples_sheet Either string or data.frame.
#'                      string: path to the TSV file
#'                      data.frame: the sample sheet
#' @param column String. Column containing the sample IDs
#' @param verbose Boolean. Whether to print the sample sheet and the IDs.
#'
#' @return Vector of strings. IDs of the samples.
#'
#' @examples
#'   samples <- get_sample_IDs("~/project/samples.tsv", "PatientID")
#'
get_sample_IDs <- function(samples_sheet, column="SampleID", verbose=TRUE) {
  if (is.character(samples_sheet)) {
    samples_metadata <- read_samples_sheet(samples_sheet, column, verbose)
  } else {
    if (!is.data.frame(samples_sheet))
      stop("samples_sheet must be a data.frame or a string")
    samples_metadata <- samples_sheet
  }
  samples <- subset(samples_metadata, select = c(column))
  samples <- samples[[1]]

  if (verbose)
    print(paste("Samples:", paste(samples, collapse=", ")))

  samples
}

