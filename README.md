# rmztabm2mtbls

R implementation of the [EBI-MetaboLights mztabm2mtbls](https://github.com/EBI-Metabolights/mztabm2mtbls) converter. It turns mzTab-M files (or their JSON export) into ISA-Tab files ready for MetaboLights submission.

## What it does
- Optionally calls the `jmztab-m` BioContainer (via Docker/Podman) to turn an mzTab-M text file into JSON.
- Maps mzTab-M metadata, samples, assays, and small molecule summaries into ISA-Tab investigation, sample, assay, and MAF files using bundled templates.
- Drops the ISA-Tab files into `output/` by default (`i_Investigation.txt`, sample sheet, assay sheet, and MAF).

## Quick start
```bash
# Install R deps
Rscript -e 'install.packages(c("jsonlite", "digest"))'

# Convert (uses Docker by default to create the JSON if needed)
Rscript convert.R \
  --input-file path/to/file.mzTab \
  --mtbls_accession_number MTBLS100001 \
  --output_dir output
```

If you already have the JSON form of the mzTab-M file:
```bash
Rscript convert.R --input-file path/to/file.mzTab.json
```

### Arguments
- `--input-file` (required): mzTab-M text or JSON file.
- `--output_dir`: target folder for ISA-Tab files (default `output`).
- `--mtbls_accession_number`: provisional study ID used to name files.
- `--container_engine`: container runtime for jmztab-m (default `docker`).
- `--mztab2m_json_convertor_image`: jmztab-m image (default `quay.io/biocontainers/jmztab-m:1.0.6--hdfd78af_1`).
- `--override_mztab2m_json_file`: force re-creation of the JSON (flag).
- `--mztabm_validation_level`: jmztab-m validation level (`Error` / `Warn` / `Info`).
- `--mztabm_mapping_file`: optional mapping XML for semantic validation.
- `--temp_folder`: temp working folder (default `output/temp`).

## Notes and scope
- Templates come from `inst/resources` and mirror the upstream Python project.
- This R port focuses on the main metadata, sample, assay, and small molecule summary mappings; less-common mzTab-M fields may need additional handling.
- Docker/Podman must be available if you want automatic mzTab-M→JSON conversion; otherwise provide the JSON yourself.
