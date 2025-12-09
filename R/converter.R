calculate_sha256 <- function(path) {
  if (!file.exists(path)) {
    return("")
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Please install the 'digest' package to compute file hashes.")
  }
  digest::digest(file = path, algo = "sha256")
}

load_mztab_json <- function(path) {
  data <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  replace_null_strings(data)
}

build_samples <- function(mztab) {
  tbl <- read_template_table(get_resource_path("s_MTBLS.txt"))
  tbl <- add_single_column(tbl, "Comment[mztab:metadata:sample:id]")
  tbl <- add_single_column(tbl, "Comment[mztab:metadata:sample:description]")

  metadata <- mztab$metadata
  samples <- metadata$sample %||% list()
  assays <- metadata$assay %||% list()
  assays_map <- list()
  for (assay in assays) {
    assays_map[[as.character(assay$id)]] <- assay
  }

  custom_characteristics <- character(0)
  for (sample in samples) {
    if (!is.null(sample$cell_type) && length(sample$cell_type) > 0) {
      custom_characteristics <- unique(c(custom_characteristics, "Cell type"))
    }
    if (!is.null(sample$disease) && length(sample$disease) > 0) {
      custom_characteristics <- unique(c(custom_characteristics, "Disease"))
    }
    if (!is.null(sample$custom)) {
      for (param in sample$custom) {
        if (!is.null(param$name)) {
          custom_characteristics <- unique(c(custom_characteristics, param$name))
        }
      }
    }
  }
  # study variables map to factors
  factor_names <- character(0)
  sample_factors <- list()
  if (!is.null(metadata$study_variable)) {
    for (sv in metadata$study_variable) {
      if (is.null(sv$factors)) next
      for (factor in sv$factors) {
        fname <- sanitize_scalar(factor$name)
        if (nzchar(fname)) {
          factor_names <- unique(c(factor_names, fname))
        }
        if (!is.null(sv$assay_refs)) {
          for (assay_ref in sv$assay_refs) {
            assay <- assays_map[[as.character(assay_ref)]]
            if (is.null(assay) || is.null(assay$sample_ref)) next
            sid <- as.character(assay$sample_ref)
            sample_factors[[sid]] <- sample_factors[[sid]] %||% list()
            val <- sanitize_scalar(sv$name %||% factor$value)
            sample_factors[[sid]][[fname]] <- unique(c(sample_factors[[sid]][[fname]] %||% character(0), val))
          }
        }
      }
    }
  }

  # add custom characteristic columns
  for (characteristic in custom_characteristics) {
    tbl <- add_ontology_columns(tbl, sprintf("Characteristics[%s]", characteristic))
  }
  for (fname in factor_names) {
    tbl <- add_ontology_columns(tbl, sprintf("Factor Value[%s]", fname))
  }

  tbl <- add_rows(tbl, length(samples))
  proto_idx <- find_header_index(tbl, "Protocol REF")
  if (!is.na(proto_idx) && nrow(tbl$data) > 0) {
    tbl$data[[proto_idx]] <- rep("Sample collection", nrow(tbl$data))
  }
  for (row_idx in seq_along(samples)) {
    sample <- samples[[row_idx]]
    species <- first_item(sample$species)
    tissue <- first_item(sample$tissue)
    tbl <- set_ontology_cell(tbl, "Characteristics[Organism]", row_idx, species$name %||% "", species$cv_label %||% "", convert_accession_number(species$cv_label, species$cv_accession))
    tbl <- set_ontology_cell(tbl, "Characteristics[Organism part]", row_idx, tissue$name %||% "", tissue$cv_label %||% "", convert_accession_number(tissue$cv_label, tissue$cv_accession))
    tbl <- set_cell(tbl, "Sample Name", row_idx, sample$name %||% "")
    tbl <- set_cell(tbl, "Source Name", row_idx, sample$name %||% "")
    tbl <- set_cell(tbl, "Comment[mztab:metadata:sample:id]", row_idx, sample$id %||% "")
    tbl <- set_cell(tbl, "Comment[mztab:metadata:sample:description]", row_idx, sample$description %||% "")
    if (!is.null(sample$disease) && length(sample$disease) > 0) {
      dis <- first_item(sample$disease)
      tbl <- set_ontology_cell(tbl, "Characteristics[Disease]", row_idx, dis$name %||% "", dis$cv_label %||% "", convert_accession_number(dis$cv_label, dis$cv_accession))
    }
    if (!is.null(sample$cell_type) && length(sample$cell_type) > 0) {
      ct <- first_item(sample$cell_type)
      tbl <- set_ontology_cell(tbl, "Characteristics[Cell type]", row_idx, ct$name %||% "", ct$cv_label %||% "", convert_accession_number(ct$cv_label, ct$cv_accession))
    }
    if (!is.null(sample$custom)) {
      for (param in sample$custom) {
        header <- sprintf("Characteristics[%s]", param$name %||% "")
        tbl <- set_ontology_cell(tbl, header, row_idx, param$value %||% param$name %||% "", param$cv_label %||% "", convert_accession_number(param$cv_label, param$cv_accession))
      }
    }
    sid <- as.character(sample$id %||% "")
    if (!is.null(sample_factors[[sid]])) {
      for (fname in names(sample_factors[[sid]])) {
        header <- sprintf("Factor Value[%s]", fname)
        tbl <- set_cell(tbl, header, row_idx, paste(sample_factors[[sid]][[fname]], collapse = ";"))
      }
    }
  }
  tbl
}

build_assay <- function(mztab, study_id) {
  tbl <- read_template_table(get_resource_path("a_template_MS-phase_metabolite_profiling.txt"), has_default_row = TRUE)
  new_column_index <- 1
  for (header in c(
    "Comment[mztab:metadata:assay:id]",
    "Comment[mztab:metadata:sample:id]",
    "Comment[mztab:metadata:ms_run:id]",
    "Comment[mztab:metadata:ms_run:name]",
    "Comment[mztab:metadata:ms_run:instrument:id]"
  )) {
    tbl <- add_single_column(tbl, header, index = new_column_index)
    new_column_index <- new_column_index + 1
  }

  metadata <- mztab$metadata
  assays <- metadata$assay %||% list()
  samples_map <- list()
  for (s in metadata$sample %||% list()) {
    samples_map[[as.character(s$id)]] <- s
  }
  ms_run_map <- list()
  for (ms in metadata$ms_run %||% list()) {
    ms_run_map[[as.character(ms$id)]] <- ms
  }
  instrument_map <- list()
  for (inst in metadata$instrument %||% list()) {
    instrument_map[[as.character(inst$id)]] <- inst
  }

  # Check whether we need extra columns
  add_checksum <- FALSE
  add_checksum_type <- FALSE
  add_raw_format <- FALSE
  add_native_id <- FALSE
  for (ms in metadata$ms_run %||% list()) {
    if (!is.null(ms$hash) && nzchar(ms$hash)) add_checksum <- TRUE
    if (!is.null(ms$hash_method) && !is.null(ms$hash_method$name)) add_checksum_type <- TRUE
    if (!is.null(ms$format) && !is.null(ms$format$name)) add_raw_format <- TRUE
    if (!is.null(ms$id_format) && !is.null(ms$id_format$name)) add_native_id <- TRUE
  }
  if (add_checksum) {
    tbl <- add_single_column(tbl, "Parameter Value[Data file checksum]")
  }
  if (add_checksum_type) {
    tbl <- add_ontology_columns(tbl, "Parameter Value[Data file checksum type]")
  }
  if (add_native_id) {
    tbl <- add_ontology_columns(tbl, "Parameter Value[Native spectrum identifier format]")
  }
  if (add_raw_format) {
    tbl <- add_ontology_columns(tbl, "Parameter Value[Raw data file format]")
  }

  tbl <- add_rows(tbl, length(assays))
  maf_name <- sprintf("m_%s_metabolite_profiling_v2_maf.tsv", study_id)

  for (row_idx in seq_along(assays)) {
    assay <- assays[[row_idx]]
    sample <- samples_map[[as.character(assay$sample_ref)]] %||% list()
    ms_run <- if (!is.null(assay$ms_run_ref) && length(assay$ms_run_ref) > 0) {
      ms_run_map[[as.character(assay$ms_run_ref[[1]])]]
    } else {
      NULL
    }
    tbl <- set_cell(tbl, "Sample Name", row_idx, sample$name %||% "")
    tbl <- set_cell(tbl, "Extract Name", row_idx, sample$name %||% "")
    tbl <- set_cell(tbl, "Labeled Extract Name", row_idx, sample$name %||% "")
    tbl <- set_cell(tbl, "MS Assay Name", row_idx, assay$name %||% "")
    tbl <- set_cell(tbl, "Metabolite Assignment File", row_idx, maf_name)
    tbl <- set_cell(tbl, "Comment[mztab:metadata:assay:id]", row_idx, assay$id %||% "")
    tbl <- set_cell(tbl, "Comment[mztab:metadata:sample:id]", row_idx, sample$id %||% "")
    if (!is.null(ms_run)) {
      tbl <- set_cell(tbl, "Raw Spectral Data File", row_idx, basename(ms_run$location %||% ""))
      tbl <- set_cell(tbl, "Comment[mztab:metadata:ms_run:id]", row_idx, ms_run$id %||% "")
      tbl <- set_cell(tbl, "Comment[mztab:metadata:ms_run:name]", row_idx, ms_run$name %||% "")
      if (!is.null(ms_run$instrument_ref)) {
        tbl <- set_cell(tbl, "Comment[mztab:metadata:ms_run:instrument:id]", row_idx, ms_run$instrument_ref %||% "")
        inst <- instrument_map[[as.character(ms_run$instrument_ref)]]
        if (!is.null(inst)) {
          tbl <- set_ontology_cell(tbl, "Parameter Value[Instrument]", row_idx, inst$name$name %||% "", inst$name$cv_label %||% "", convert_accession_number(inst$name$cv_label, inst$name$cv_accession))
          tbl <- set_ontology_cell(tbl, "Parameter Value[Ion source]", row_idx, inst$source$name %||% "", inst$source$cv_label %||% "", convert_accession_number(inst$source$cv_label, inst$source$cv_accession))
          if (!is.null(inst$analyzer) && length(inst$analyzer) > 0) {
            an <- first_item(inst$analyzer)
            tbl <- set_ontology_cell(tbl, "Parameter Value[Mass analyzer]", row_idx, an$name %||% "", an$cv_label %||% "", convert_accession_number(an$cv_label, an$cv_accession))
          }
        }
      }
      if (!is.null(ms_run$scan_polarity) && length(ms_run$scan_polarity) > 0) {
        pol <- first_item(ms_run$scan_polarity)
        tbl <- set_ontology_cell(tbl, "Parameter Value[Scan polarity]", row_idx, pol$name %||% "", pol$cv_label %||% "", convert_accession_number(pol$cv_label, pol$cv_accession))
      }
      if (add_checksum) {
        tbl <- set_cell(tbl, "Parameter Value[Data file checksum]", row_idx, ms_run$hash %||% "")
      }
      if (add_checksum_type) {
        h <- ms_run$hash_method %||% list()
        tbl <- set_ontology_cell(tbl, "Parameter Value[Data file checksum type]", row_idx, h$name %||% "", h$cv_label %||% "", convert_accession_number(h$cv_label, h$cv_accession))
      }
      if (add_native_id) {
        idf <- ms_run$id_format %||% list()
        tbl <- set_ontology_cell(tbl, "Parameter Value[Native spectrum identifier format]", row_idx, idf$name %||% "", idf$cv_label %||% "", convert_accession_number(idf$cv_label, idf$cv_accession))
      }
      if (add_raw_format) {
        fmt <- ms_run$format %||% list()
        tbl <- set_ontology_cell(tbl, "Parameter Value[Raw data file format]", row_idx, fmt$name %||% "", fmt$cv_label %||% "", convert_accession_number(fmt$cv_label, fmt$cv_accession))
      }
    }
  }
  tbl
}

build_maf <- function(mztab, assays) {
  tbl <- read_template_table(get_resource_path("m_MTBLS_metabolite_profiling_v2_maf.tsv"))
  tbl <- add_single_column(tbl, "Comment[mztab:summary:sml_id]", index = 1)
  for (idx in seq_along(assays)) {
    tbl <- add_single_column(tbl, sanitize_scalar(assays[[idx]]$name %||% paste0("assay_", idx)))
  }
  summaries <- mztab$smallMoleculeSummary %||% list()
  features <- mztab$smallMoleculeFeature %||% list()
  feature_map <- list()
  for (feat in features) {
    feature_map[[as.character(feat$smf_id)]] <- feat
  }
  tbl <- add_rows(tbl, length(summaries))
  databases <- mztab$metadata$database %||% list()
  db_prefixes <- vapply(databases, function(x) sanitize_scalar(x$prefix), character(1))
  for (row_idx in seq_along(summaries)) {
    sms <- summaries[[row_idx]]
    setter <- function(header, value) {
      tbl <<- set_cell(tbl, header, row_idx, value)
    }
    setter("database_identifier", paste(sanitize_value(sms$database_identifier %||% list()), collapse = "|"))
    setter("chemical_formula", paste(sanitize_value(sms$chemical_formula %||% list()), collapse = "|"))
    setter("smiles", paste(sanitize_value(sms$smiles %||% list()), collapse = "|"))
    setter("inchi", paste(sanitize_value(sms$inchi %||% list()), collapse = "|"))
    setter("metabolite_identification", paste(sanitize_value(sms$chemical_name %||% list()), collapse = "|"))
    setter("reliability", sanitize_scalar(sms$best_id_confidence_value %||% sms$reliability %||% ""))
    setter("Comment[mztab:summary:sml_id]", sms$sml_id %||% "")
    if (!is.null(sms$smf_id_refs)) {
      rt_vals <- c()
      for (ref in sms$smf_id_refs) {
        feat <- feature_map[[as.character(ref)]]
        if (!is.null(feat$exp_mass_to_charge)) {
          rt_vals <- c(rt_vals, feat$exp_mass_to_charge)
        }
      }
      setter("retention_time", paste(sanitize_value(rt_vals), collapse = "|"))
    }
    if (!is.null(sms$theoretical_neutral_mass)) {
      setter("mass_to_charge", paste(sanitize_value(sms$theoretical_neutral_mass), collapse = "|"))
    }
    setter("database", paste(db_prefixes[db_prefixes != ""], collapse = ";"))
    setter("database_version", paste(db_prefixes[db_prefixes != ""], collapse = ";"))
    if (!is.null(sms$abundance_assay)) {
      start_col <- 2 # after comment column
      for (i in seq_along(sms$abundance_assay)) {
        col_idx <- start_col + i - 1
        if (col_idx <= length(tbl$headers)) {
          tbl$data[[col_idx]][row_idx] <- sanitize_scalar(sms$abundance_assay[[i]])
        }
      }
    }
  }
  tbl
}

convert_mztabm <- function(input_file,
                           output_dir = "output",
                           mtbls_accession_number = "MTBLS1000000",
                           container_engine = "docker",
                           mztab2m_json_convertor_image = "quay.io/biocontainers/jmztab-m:1.0.6--hdfd78af_1",
                           override_mztab2m_json_file = FALSE,
                           mztabm_validation_level = "Error",
                           mztabm_mapping_file = NULL,
                           temp_folder = NULL) {
  if (is.null(temp_folder)) {
    temp_folder <- file.path(output_dir, "temp")
  }
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)

  if (missing(input_file) || is.null(input_file) || !file.exists(input_file)) {
    stop("Please provide an existing mzTab-M file path via input_file.")
  }
  message("'", input_file, "' will be converted to MetaboLights ISA-TAB metadata files.")
  source_hash <- calculate_sha256(input_file)
  input_json_file <- input_file
  ext <- tolower(tools::file_ext(input_file))
  if (ext != "json") {
    input_json_file <- file.path(temp_folder, paste0(basename(input_file), ".json"))
    if (!override_mztab2m_json_file && file.exists(input_json_file)) {
      message("Reusing existing mzTab-M JSON file: ", input_json_file)
    } else {
      file.copy(input_file, file.path(temp_folder, basename(input_file)), overwrite = TRUE)
      success <- run_jmztabm_docker(
        container_engine = container_engine,
        mztab2m_json_convertor_image = mztab2m_json_convertor_image,
        dirname = temp_folder,
        filename = basename(input_file),
        mztabm_validation_level = mztabm_validation_level,
        mztabm_mapping_file = mztabm_mapping_file
      )
      if (!success) {
        stop("mzTab-M to JSON conversion failed.")
      }
    }
  }
  if (!file.exists(input_json_file)) {
    stop("JSON input file not found at ", input_json_file)
  }
  mztab <- load_mztab_json(input_json_file)

  inv <- build_investigation(mztab, mtbls_accession_number, input_file, source_hash)
  samples_tbl <- build_samples(mztab)
  assay_tbl <- build_assay(mztab, mtbls_accession_number)
  maf_tbl <- build_maf(mztab, mztab$metadata$assay %||% list())

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_investigation(inv, file.path(output_dir, "i_Investigation.txt"))
  write_table(samples_tbl, file.path(output_dir, sprintf("s_%s.txt", mtbls_accession_number)))
  write_table(assay_tbl, file.path(output_dir, sprintf("a_%s_metabolite_profiling.txt", mtbls_accession_number)))
  write_table(maf_tbl, file.path(output_dir, sprintf("m_%s_metabolite_profiling_v2_maf.tsv", mtbls_accession_number)))
  invisible(list(
    investigation = file.path(output_dir, "i_Investigation.txt"),
    samples = file.path(output_dir, sprintf("s_%s.txt", mtbls_accession_number)),
    assay = file.path(output_dir, sprintf("a_%s_metabolite_profiling.txt", mtbls_accession_number)),
    maf = file.path(output_dir, sprintf("m_%s_metabolite_profiling_v2_maf.tsv", mtbls_accession_number))
  ))
}
