read_investigation_template <- function() {
  path <- get_resource_path("i_Investigation.txt")
  read.delim(path,
    header = FALSE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = ""
  )
}

ensure_columns <- function(df, needed) {
  if (ncol(df) >= needed) {
    return(df)
  }
  extra <- needed - ncol(df)
  for (i in seq_len(extra)) {
    df[[ncol(df) + 1]] <- ""
  }
  df
}

set_investigation_field <- function(df, field, values) {
  idx <- which(df[[1]] == field)
  if (length(idx) == 0) {
    new_row <- c(field, rep("", max(length(values), 1)))
    df <- rbind(df, new_row)
    idx <- nrow(df)
  }
  df <- ensure_columns(df, length(values) + 1)
  df[idx, 2:(length(values) + 1)] <- sanitize_value(values)
  df
}

append_comment_row <- function(df, name, values) {
  new_row <- c(name, sanitize_value(values))
  df <- ensure_columns(df, length(new_row))
  df <- rbind(df, new_row)
  df
}

add_cv_to_investigation <- function(df, cv_list) {
  if (is.null(cv_list) || length(cv_list) == 0) {
    return(df)
  }
  fields <- c(
    "Term Source Name",
    "Term Source File",
    "Term Source Version",
    "Term Source Description",
    "Comment[mztab:metadata:cv:id]"
  )
  current_cols <- ncol(df)
  for (field in fields) {
    if (!field %in% df[[1]]) {
      df <- rbind(df, c(field, rep("", current_cols - 1)))
    }
  }
  for (idx in seq_along(cv_list)) {
    cv <- cv_list[[idx]]
    label <- sanitize_scalar(cv$label)
    uri <- sanitize_scalar(cv$uri)
    version <- sanitize_scalar(cv$version)
    desc <- sanitize_scalar(cv$full_name)
    id <- sanitize_scalar(cv$id)
    df <- ensure_columns(df, ncol(df) + 1)
    col <- ncol(df)
    for (field in fields) {
      row_idx <- which(df[[1]] == field)
      val <- switch(field,
        "Term Source Name" = label,
        "Term Source File" = uri,
        "Term Source Version" = version,
        "Term Source Description" = desc,
        "Comment[mztab:metadata:cv:id]" = id,
        ""
      )
      df[row_idx, col] <- val
    }
  }
  df
}

add_contacts_to_investigation <- function(df, contacts) {
  if (is.null(contacts) || length(contacts) == 0) {
    return(df)
  }
  rows <- c(
    "Study Person Last Name",
    "Study Person First Name",
    "Study Person Mid Initials",
    "Study Person Email",
    "Study Person Affiliation",
    "Study Person Roles"
  )
  for (row in rows) {
    if (!row %in% df[[1]]) {
      df <- rbind(df, c(row, rep("", ncol(df) - 1)))
    }
  }
  col_offset <- max(2, ncol(df) + 1)
  df <- ensure_columns(df, col_offset - 1)
  for (contact in contacts) {
    name <- sanitize_scalar(contact$name)
    parts <- strsplit(name, " +")[[1]]
    first <- if (length(parts) > 0) parts[[1]] else ""
    last <- if (length(parts) > 1) parts[[length(parts)]] else ""
    mid <- if (length(parts) > 2) paste(parts[2:(length(parts) - 1)], collapse = " ") else ""
    df <- ensure_columns(df, ncol(df) + 1)
    col_idx <- ncol(df)
    df[df[[1]] == "Study Person First Name", col_idx] <- first
    df[df[[1]] == "Study Person Last Name", col_idx] <- last
    df[df[[1]] == "Study Person Mid Initials", col_idx] <- mid
    df[df[[1]] == "Study Person Email", col_idx] <- sanitize_scalar(contact$email)
    df[df[[1]] == "Study Person Affiliation", col_idx] <- sanitize_scalar(contact$affiliation)
    role <- if (col_idx == 2) "Principal Investigator" else "Author"
    df[df[[1]] == "Study Person Roles", col_idx] <- role
  }
  df
}

add_publications_to_investigation <- function(df, publications, mztab) {
  if (is.null(publications) || length(publications) == 0) {
    return(df)
  }
  dois <- character(0)
  pubmeds <- character(0)
  uris <- character(0)
  ids <- character(0)
  for (pub in publications) {
    if (!is.null(pub$publicationItems)) {
      for (item in pub$publicationItems) {
        if (identical(item$type, "doi")) {
          dois <- c(dois, sanitize_scalar(item$accession))
        } else if (identical(item$type, "pubmed")) {
          pubmeds <- c(pubmeds, sanitize_scalar(item$accession))
        } else if (identical(item$type, "uri")) {
          uris <- c(uris, sanitize_scalar(item$accession))
        }
      }
    }
    ids <- c(ids, sanitize_scalar(pub$id))
  }
  df <- set_investigation_field(df, "Study Publication DOI", paste(dois, collapse = ";"))
  df <- set_investigation_field(df, "Study PubMed ID", paste(pubmeds, collapse = ";"))
  df <- set_investigation_field(df, "Comment[mztab:metadata:publication:uri]", paste(uris, collapse = ";"))
  df <- set_investigation_field(df, "Comment[mztab:metadata:publication:id]", paste(ids, collapse = ";"))
  if (!is.null(mztab$metadata$title)) {
    df <- set_investigation_field(df, "Study Publication Title", sanitize_scalar(mztab$metadata$title))
  }
  if (!is.null(mztab$metadata$contact)) {
    authors <- vapply(mztab$metadata$contact, function(x) sanitize_scalar(x$name), character(1))
    df <- set_investigation_field(df, "Study Publication Author List", paste(authors, collapse = ", "))
  }
  df
}

build_investigation <- function(mztab, study_id, source_file, source_hash) {
  inv <- read_investigation_template()
  inv <- set_investigation_field(inv, "Investigation Identifier", study_id)
  inv <- set_investigation_field(inv, "Study Identifier", study_id)
  inv <- set_investigation_field(inv, "Study File Name", sprintf("s_%s.txt", study_id))
  inv <- set_investigation_field(inv, "Study Assay File Name", sprintf("a_%s_metabolite_profiling.txt", study_id))
  inv <- set_investigation_field(inv, "Study Title", sanitize_scalar(mztab$metadata$title))
  inv <- set_investigation_field(inv, "Study Description", sanitize_scalar(mztab$metadata$description))
  today <- format(Sys.Date(), "%Y-%m-%d")
  inv <- set_investigation_field(inv, "Study Submission Date", today)
  inv <- set_investigation_field(inv, "Study Public Release Date", today)
  inv <- set_investigation_field(inv, "Investigation Submission Date", today)
  inv <- set_investigation_field(inv, "Investigation Public Release Date", today)
  inv <- set_investigation_field(inv, "Comment[mztab:metadata:mzTab_version]", sanitize_scalar(mztab$metadata$mzTab_version))
  inv <- set_investigation_field(inv, "Comment[mztab:metadata:mzTab_ID]", sanitize_scalar(mztab$metadata$mzTab_ID))
  if (!is.null(mztab$metadata$uri)) {
    uris <- vapply(mztab$metadata$uri, function(x) sanitize_scalar(x$value %||% x), character(1))
    inv <- set_investigation_field(inv, "Comment[mztab:metadata:uri]", paste(uris, collapse = ";"))
  }
  if (!is.null(mztab$metadata$external_study_uri)) {
    uris <- vapply(mztab$metadata$external_study_uri, function(x) sanitize_scalar(x$value %||% x), character(1))
    inv <- set_investigation_field(inv, "Comment[mztab:metadata:external_study_uri]", paste(uris, collapse = ";"))
  }
  inv <- append_comment_row(inv, "Comment[mztab:source_file:location]", sanitize_scalar(source_file))
  inv <- append_comment_row(inv, "Comment[mztab:source_file:hash:sha256]", sanitize_scalar(source_hash))
  inv <- add_cv_to_investigation(inv, mztab$metadata$cv)
  inv <- add_contacts_to_investigation(inv, mztab$metadata$contact)
  inv <- add_publications_to_investigation(inv, mztab$metadata$publication, mztab)
  inv
}

write_investigation <- function(inv, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(inv, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
}
