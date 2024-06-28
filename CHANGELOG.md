## SINCLAIR development version

- Overhaul the CLI to use python rather than bash, which introduces breaking changes (@kelly-sovacool #61).
  - Create a script (`bin/sinclair`) to provide an interface to the CLI that works out-of-the-box without the need to install the python package with pip. (#80, @kelly-sovacool)
- Bug fixes:
  - Fix biowulf module syntax in `conf/modules.config`. (@epehrsson, #81)
  - Fix filtering thresholds and use filtered object for downstream steps in `SEURAT_PROCESS`. (@epehrsson, #81)
  - Fix seurat object and group assignment in `SEURAT_MERGE`. (@epehrsson, #81)

## SINCLAIR 0.2.0

TODO fill in prior changelog entries from github release notes

## SINCLAIR 0.1.1

## SINCLAIR 0.1.0
