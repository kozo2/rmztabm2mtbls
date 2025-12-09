run_jmztabm_docker <- function(container_engine = "docker",
                               mztab2m_json_convertor_image = "quay.io/biocontainers/jmztab-m:1.0.6--hdfd78af_1",
                               dirname = ".",
                               filename = NULL,
                               mztabm_validation_level = "Error",
                               mztabm_mapping_file = NULL) {
  if (is.null(filename)) {
    stop("filename is required to run jmztab-m")
  }

  local_cmd <- c(container_engine, "run", "--rm")
  mounts <- c("-v", sprintf("%s:/home/data", normalizePath(dirname)))
  mapping_args <- character(0)
  if (!is.null(mztabm_mapping_file)) {
    mapping_abs <- normalizePath(mztabm_mapping_file)
    mapping_name <- basename(mapping_abs)
    mounts <- c(mounts, "-v", sprintf("%s:/home/configuration/%s", mapping_abs, mapping_name))
    mapping_args <- c("-s", sprintf("/home/configuration/%s", mapping_name))
  }
  jmztab_cmd <- c(
    "--workdir=/home/data",
    mztab2m_json_convertor_image,
    "jmztab-m",
    "-c",
    sprintf("/home/data/%s", filename),
    "--toJson",
    "-o",
    sprintf("/home/data/%s.validation.txt", filename),
    "-l",
    mztabm_validation_level
  )
  jmztab_cmd <- c(jmztab_cmd, mapping_args)
  full_cmd <- c(local_cmd, mounts, jmztab_cmd)
  message("Running jmztab-m: ", paste(full_cmd, collapse = " "))
  res <- system2(full_cmd[1], args = full_cmd[-1], stdout = TRUE, stderr = TRUE)
  attr(res, "status") <- attr(res, "status") %||% 0
  status <- attr(res, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    warning("jmztab-m docker run failed with status ", status, call. = FALSE)
    return(FALSE)
  }
  TRUE
}
