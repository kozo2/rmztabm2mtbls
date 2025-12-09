sanitize_value <- function(value) {
  if (is.null(value)) {
    return("")
  }
  if (is.list(value)) {
    return(lapply(value, sanitize_value))
  }
  gsub("[\r\n\t]", " ", trimws(as.character(value)))
}

sanitize_scalar <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("")
  }
  gsub("[\r\n\t]", " ", trimws(as.character(value[1])))
}

"%||%" <- function(a, b) {
  if (is.null(a)) b else a
}

replace_null_strings <- function(x) {
  if (is.list(x)) {
    for (idx in seq_along(x)) {
      x[[idx]] <- replace_null_strings(x[[idx]])
    }
    return(x)
  }
  if (is.character(x) && length(x) == 1 && identical(x, "null")) {
    return(NULL)
  }
  x
}

convert_accession_number <- function(cv_label, cv_accession) {
  acc <- sanitize_scalar(cv_accession)
  if (acc == "") {
    return("")
  }
  parts <- strsplit(acc, ":", fixed = TRUE)[[1]]
  accession_number <- if (length(parts) > 1) parts[[2]] else parts[[1]]
  sprintf("http://purl.obolibrary.org/obo/%s_%s", sanitize_scalar(cv_label), accession_number)
}

get_resource_path <- function(filename) {
  sys_path <- system.file("resources", filename, package = "rmztabm2mtbls")
  if (nzchar(sys_path)) {
    return(sys_path)
  }
  file.path("inst", "resources", filename)
}

first_item <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(list())
  }
  x[[1]]
}

read_template_table <- function(path, has_default_row = FALSE) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    stop(sprintf("Template %s is empty", path))
  }
  header <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
  default_row <- NULL
  if (has_default_row && length(lines) > 1) {
    default_row <- strsplit(lines[[2]], "\t", fixed = TRUE)[[1]]
  }
  make_table(header, default_row = default_row)
}

make_table <- function(headers, default_row = NULL) {
  col_count <- length(headers)
  data <- data.frame(matrix(character(), nrow = 0, ncol = col_count), stringsAsFactors = FALSE)
  colnames(data) <- paste0("col", seq_len(col_count))
  list(
    headers = headers,
    data = data,
    default_row = if (!is.null(default_row)) default_row else rep("", col_count)
  )
}

sync_colnames <- function(tbl) {
  colnames(tbl$data) <- paste0("col", seq_len(ncol(tbl$data)))
  tbl
}

add_rows <- function(tbl, n) {
  if (n <= 0) {
    return(tbl)
  }
  base_row <- tbl$default_row
  if (length(base_row) < length(tbl$headers)) {
    base_row <- c(base_row, rep("", length(tbl$headers) - length(base_row)))
  } else if (length(base_row) > length(tbl$headers)) {
    base_row <- base_row[seq_along(tbl$headers)]
  }
  new_data <- matrix(rep(base_row, n), nrow = n, byrow = TRUE)
  colnames(new_data) <- colnames(tbl$data)
  tbl$data <- rbind(tbl$data, as.data.frame(new_data, stringsAsFactors = FALSE))
  tbl
}

insert_column <- function(tbl, header_label, index = NULL, fill = "") {
  if (is.null(index) || index < 1 || index > length(tbl$headers) + 1) {
    index <- length(tbl$headers) + 1
  }
  tbl$headers <- append(tbl$headers, header_label, after = index - 1)

  fill_default <- fill
  tbl$default_row <- append(tbl$default_row, fill_default, after = index - 1)

  if (ncol(tbl$data) == 0) {
    tbl$data <- data.frame(matrix(character(), nrow = 0, ncol = length(tbl$headers)), stringsAsFactors = FALSE)
  } else {
    left <- if (index > 1) tbl$data[, seq_len(index - 1), drop = FALSE] else NULL
    right <- if (index <= ncol(tbl$data)) tbl$data[, seq(index, ncol(tbl$data)), drop = FALSE] else NULL
    new_col <- rep(fill, nrow(tbl$data))
    tbl$data <- cbind(left, new_col, right)
  }
  tbl <- sync_colnames(tbl)
  tbl
}

add_single_column <- function(tbl, header_label, index = NULL, fill = "") {
  insert_column(tbl, header_label, index = index, fill = fill)
}

add_ontology_columns <- function(tbl, header_label, index = NULL) {
  tbl <- insert_column(tbl, header_label, index = index, fill = "")
  tbl <- insert_column(tbl, "Term Source REF", index = if (is.null(index)) NULL else index + 1, fill = "")
  tbl <- insert_column(tbl, "Term Accession Number", index = if (is.null(index)) NULL else index + 2, fill = "")
  tbl
}

find_header_index <- function(tbl, header_label) {
  match(header_label, tbl$headers)
}

set_cell <- function(tbl, header_label, row, value) {
  idx <- find_header_index(tbl, header_label)
  if (is.na(idx)) {
    return(tbl)
  }
  if (row > nrow(tbl$data)) {
    stop("Row index out of range when setting value")
  }
  tbl$data[[idx]][row] <- sanitize_scalar(value)
  tbl
}

set_ontology_cell <- function(tbl, header_label, row, value, source_ref = "", accession = "") {
  idx <- find_header_index(tbl, header_label)
  if (is.na(idx)) {
    return(tbl)
  }
  tbl$data[[idx]][row] <- sanitize_scalar(value)
  if (!is.na(idx + 1) && (idx + 1) <= length(tbl$headers)) {
    tbl$data[[idx + 1]][row] <- sanitize_scalar(source_ref)
  }
  if (!is.na(idx + 2) && (idx + 2) <= length(tbl$headers)) {
    tbl$data[[idx + 2]][row] <- sanitize_scalar(accession)
  }
  tbl
}

write_table <- function(tbl, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  header_line <- paste(tbl$headers, collapse = "\t")
  body <- character(nrow(tbl$data))
  if (nrow(tbl$data) > 0) {
    for (i in seq_len(nrow(tbl$data))) {
      row_vals <- vapply(tbl$data[i, ], sanitize_scalar, character(1))
      body[[i]] <- paste(row_vals, collapse = "\t")
    }
  }
  writeLines(c(header_line, body), con = path)
}
