# CCBR/SINCLAIR pipeline parameters

SINgle CelL AnalysIs Resource

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter        | Description                                                                                                                                                                                                                                                                                                                                                                      | Type      | Default                      | Required |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ---------------------------- | -------- |
| `input`          | Path to CSV file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row.</small></details> | `string`  | assets/input_manifest.csv    | Yes      |
| `contrast`       | Path to CSV file with contrast specification.                                                                                                                                                                                                                                                                                                                                    | `string`  | assets/contrast_manifest.csv | Yes      |
| `outdir`         | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.                                                                                                                                                                                                                                                         | `string`  | output                       | Yes      |
| `run_cellranger` | Whether to run Cell Ranger.                                                                                                                                                                                                                                                                                                                                                      | `boolean` | True                         | No       |

## Main options

| Parameter            | Description                                                                                           | Type      | Default                   |
| -------------------- | ----------------------------------------------------------------------------------------------------- | --------- | ------------------------- |
| `species`            | Reference genome species to use.                                                                      | `string`  | hg19                      |
| `vars_to_regress`    | Comma-separated list of variables to regress out during data scaling (e.g., `nCount_RNA,percent_mt`). | `string`  |                           |
| `qc_filtering`       | Type of QC filtering to perform (e.g., `manual`, `auto`).                                             | `string`  | manual                    |
| `nCount_RNA_max`     | Maximum total UMI count per cell for QC filtering.                                                    | `integer` | 500000                    |
| `nCount_RNA_min`     | Minimum total UMI count per cell for QC filtering.                                                    | `integer` | 1000                      |
| `nFeature_RNA_max`   | Maximum number of detected genes per cell for QC filtering.                                           | `integer` | 5000                      |
| `nFeature_RNA_min`   | Minimum number of detected genes per cell for QC filtering.                                           | `integer` | 200                       |
| `percent_mt_max`     | Maximum percent of mitochondrial reads allowed per cell.                                              | `integer` | 10                        |
| `percent_mt_min`     | Minimum percent of mitochondrial reads allowed per cell.                                              | `integer` | 0                         |
| `run_doublet_finder` | Whether to run doublet detection using DoubletFinder (`Y` or `N`).                                    | `string`  | Y                         |
| `seurat_resolution`  | Comma-separated list of resolutions for Seurat clustering.                                            | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1 |
| `npcs`               | Number of principal components to use for dimensionality reduction.                                   | `integer` | 50                        |
| `resolution_list`    | Alias for `seurat_resolution`, if needed in other modules.                                            | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1 |
| `genome_dir`         | Path to the genome references directory. Overridden by platform-specific configs.                     | `string`  |                           |

## Institutional config options

| Parameter                    | Description                                                            | Type     | Default |
| ---------------------------- | ---------------------------------------------------------------------- | -------- | ------- |
| `config_profile_name`        | Name of the configuration profile (e.g., the environment or platform). | `string` |         |
| `config_profile_description` | Short description of the profile's purpose or context.                 | `string` |         |
| `config_profile_contact`     | Contact information for support related to this profile.               | `string` |         |
| `config_profile_url`         | URL for documentation or more information about the profile.           | `string` |         |

## Generic options

| Parameter          | Description                                                               | Type      | Default                        | Required |
| ------------------ | ------------------------------------------------------------------------- | --------- | ------------------------------ | -------- |
| `publish_dir_mode` | How output files are published (e.g., `copy`, `link`, or `symlink`).      | `string`  | link                           | True     |
| `tracedir`         | Directory to store pipeline execution metadata and logs.                  | `string`  | ${params.outdir}/pipeline_info | True     |
| `max_memory`       | Maximum memory that can be allocated per process.                         | `string`  | 128.GB                         | True     |
| `max_cpus`         | Maximum number of CPUs that can be allocated per process.                 | `integer` | 48                             | True     |
| `max_time`         | Maximum execution time allowed per process (e.g., `240.h` for 240 hours). | `string`  | 240.h                          | True     |

## Containers

Docker/Singularity containers to use for processes. Must be available in dockerhub

| Parameter         | Description                                                        | Type     | Default                             | Required | Hidden |
| ----------------- | ------------------------------------------------------------------ | -------- | ----------------------------------- | -------- | ------ |
| `base_container`  | Docker/Singularity image containing base environment dependencies. | `string` | nciccbr/ccbr_ubuntu_base_20.04:v6.1 |          | True   |
| `baser_container` | Container image used for running the BASER tool or module.         | `string` | nciccbr/sinclair_baser:0.1.0        |          | True   |

## Hidden options

| Parameter               | Description                                                                              | Type     | Default                                                 |
| ----------------------- | ---------------------------------------------------------------------------------------- | -------- | ------------------------------------------------------- |
| `Rlib_dir`              | Path to the directory containing the required R libraries.                               | `string` | /data/CCBR_Pipeliner/db/PipeDB/Rlibrary_4.3_scRNA_RHEL8 |
| `conda_path`            | Path to the Conda environment used for Python-based analysis.                            | `string` | /data/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env      |
| `python_path`           | Path to the Python binary within the specified Conda environment.                        | `string` | /data/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env/bin  |
| `Rpkg`                  | Path to the R package configuration file listing required packages.                      | `string` | ${projectDir}/conf/Rpack.config                         |
| `script_functions`      | R script file containing commonly used scRNA-seq helper functions.                       | `string` | ${projectDir}/bin/scRNA_functions.R                     |
| `script_preprocess`     | R Markdown script for preprocessing single-cell data with Seurat.                        | `string` | ${projectDir}/bin/seurat_preprocess.Rmd                 |
| `script_merge`          | R Markdown script for merging Seurat objects from multiple samples.                      | `string` | ${projectDir}/bin/seurat_merge.Rmd                      |
| `script_bc_harmony`     | R Markdown script for batch correction using Harmony.                                    | `string` | ${projectDir}/bin/batch_correction_harmony.Rmd          |
| `script_bc_rpca`        | R Markdown script for batch correction using Seurat’s reciprocal PCA (RPCA).             | `string` | ${projectDir}/bin/batch_correction_rpca.Rmd             |
| `script_bc_cca`         | R Markdown script for batch correction using canonical correlation analysis (CCA).       | `string` | ${projectDir}/bin/batch_correction_cca.Rmd              |
| `script_liger`          | R Markdown script for batch correction using LIGER integration.                          | `string` | ${projectDir}/bin/batch_correction_liger.Rmd            |
| `script_bc_integration` | R Markdown script for integrating batch correction results into a final unified dataset. | `string` | ${projectDir}/bin/batch_correction_integration.Rmd      |

<!-- this doc is generated by: nf-core pipelines schema docs -->
