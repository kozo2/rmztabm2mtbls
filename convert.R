#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("rmztabm2mtbls", quietly = TRUE)) {
    library(rmztabm2mtbls)
  } else {
    message("Loading rmztabm2mtbls functions from local source files.")
    r_files <- list.files("R", full.names = TRUE, pattern = "\\\\.R$")
    for (f in r_files) {
      source(f)
    }
  }
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    input_file = NULL,
    output_dir = "output",
    mtbls_accession_number = "MTBLS1000000",
    container_engine = "docker",
    mztab2m_json_convertor_image = "quay.io/biocontainers/jmztab-m:1.0.6--hdfd78af_1",
    override_mztab2m_json_file = FALSE,
    mztabm_validation_level = "Error",
    mztabm_mapping_file = NULL,
    temp_folder = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unknown argument: ", key)
    }
    key_clean <- substring(key, 3)
    if (key_clean %in% c("override_mztab2m_json_file")) {
      out[[key_clean]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    value <- args[[i + 1]]
    out[[key_clean]] <- value
    i <- i + 2
  }
  if (is.null(out$input_file)) {
    stop("Please provide --input-file <path to mzTab-M or JSON file>.")
  }
  out
}

main <- function() {
  opts <- parse_args()
  convert_mztabm(
    input_file = opts$input_file,
    output_dir = opts$output_dir,
    mtbls_accession_number = opts$mtbls_accession_number,
    container_engine = opts$container_engine,
    mztab2m_json_convertor_image = opts$mztab2m_json_convertor_image,
    override_mztab2m_json_file = opts$override_mztab2m_json_file,
    mztabm_validation_level = opts$mztabm_validation_level,
    mztabm_mapping_file = opts$mztabm_mapping_file,
    temp_folder = opts$temp_folder
  )
}

if (identical(environment(), globalenv())) {
  main()
}
