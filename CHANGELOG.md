## SINCLAIR development version

- Overhaul the CLI to use python rather than bash, which introduces breaking changes (#61, @kelly-sovacool).
  - Create a script (`bin/sinclair`) to provide an interface to the CLI that works out-of-the-box without the need to install the python package with pip. (#80, @kelly-sovacool)
- Bug fixes:
  - Fix biowulf module syntax in `conf/modules.config`. (#81, @epehrsson)
  - Fix filtering thresholds and use filtered object for downstream steps in `SEURAT_PROCESS`. (#81, @epehrsson)
  - Fix seurat object and group assignment in `SEURAT_MERGE`. (#81, @epehrsson)
  - Correctly map h5 files from cellranger with their input fastq files during preprocessing. (#82, @kelly-sovacool)
  - Fix the $SLURM_JOB_ID variable name for biowulf. (#82, @kelly-sovacool)
  - Fix file paths for test dataset. (#82, @kelly-sovacool)
  - Add the tex module for preprocess and merge, which produce PDF files from R Markdown. (#82, @kelly-sovacool)

## SINCLAIR 0.2.0

TODO fill in prior changelog entries from github release notes (see `dev` branch)

## SINCLAIR 0.1.1

## SINCLAIR 0.1.0
