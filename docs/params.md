# CCBR/SINCLAIR pipeline parameters

SINgle CelL AnalysIs Resource

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter  | Description                                                                                                                                                                                                                                                                                                                                                                                  | Type     | Default                                           | Required | Hidden |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------------------------------------------------- | -------- | ------ |
| `input`    | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details> | `string` | ${launchDir}/assets/input_manifest_cellranger.csv | True     |        |
| `contrast` |                                                                                                                                                                                                                                                                                                                                                                                              | `string` | ${launchDir}/assets/contrast_manifest.csv         |          |        |
| `outdir`   | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.                                                                                                                                                                                                                                                                     | `string` | ${launchDir}/output                               | True     |        |

## main options

Most commonly used parameters

| Parameter            | Description                                                                                        | Type      | Default                        | Required | Hidden |
| -------------------- | -------------------------------------------------------------------------------------------------- | --------- | ------------------------------ | -------- | ------ |
| `species`            |                                                                                                    | `string`  | hg19                           |          |        |
| `run_cellranger`     | whether to run cellranger                                                                          | `string`  | Y                              |          |        |
| `vars_to_regress`    | comma-separated list of variables to regress                                                       | `string`  | percent.mt,nFeature_RNA        |          |        |
| `qc_filtering`       |                                                                                                    | `string`  | manual                         |          |        |
| `nCount_RNA_max`     |                                                                                                    | `integer` | 500000                         |          |        |
| `nCount_RNA_min`     |                                                                                                    | `integer` | 1000                           |          |        |
| `nFeature_RNA_max`   |                                                                                                    | `integer` | 5000                           |          |        |
| `nFeature_RNA_min`   |                                                                                                    | `integer` | 200                            |          |        |
| `percent_mt_max`     |                                                                                                    | `integer` | 10                             |          |        |
| `percent_mt_min`     |                                                                                                    | `integer` | 0                              |          |        |
| `run_doublet_finder` |                                                                                                    | `string`  | Y                              |          |        |
| `seurat_resolution`  |                                                                                                    | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1      |          |        |
| `npcs`               | Number of principle components                                                                     | `integer` | 50                             |          |        |
| `resolution_list`    |                                                                                                    | `string`  | 0.1,0.2,0.3,0.5,0.6,0.8,1      |          |        |
| `genome_dir`         | Path to genome reference. This is set by platforms-specific config files, e.g. conf/biowulf.config | `string`  |                                |          |        |
| `tracedir`           |                                                                                                    | `string`  | ${params.outdir}/pipeline_info |          |        |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter    | Description                                                                                                                                                                                                                                                                 | Type      | Default | Required | Hidden |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `max_cpus`   | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>                                      | `integer` | 48      |          | True   |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details> | `string`  | 128.GB  |          | True   |
| `max_time`   | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>        | `string`  | 240.h   |          | True   |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter          | Description                                                                                                                                                                                                                                                                                                                                                                                                  | Type     | Default | Required | Hidden |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------- | ------- | -------- | ------ |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details> | `string` | link    |          | True   |

## Institutional config options

Parameters used to describe centralised config profiles.

| Parameter                    | Description | Type     | Default | Required | Hidden |
| ---------------------------- | ----------- | -------- | ------- | -------- | ------ |
| `config_profile_description` |             | `string` |         |          | True   |
| `config_profile_contact`     |             | `string` |         |          | True   |
| `config_profile_url`         |             | `string` |         |          | True   |

## Other parameters

| Parameter               | Description | Type     | Default                                                                                                        | Required | Hidden |
| ----------------------- | ----------- | -------- | -------------------------------------------------------------------------------------------------------------- | -------- | ------ |
| `Rlib_dir`              |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Rlibrary_4.3_scRNA_RHEL8                                           |          | True   |
| `conda_path`            |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env                                                |          | True   |
| `python_path`           |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/db/PipeDB/Conda/envs/scvi-env/bin                                            |          | True   |
| `Rpkg`                  |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/conf/Rpack.config                    |          | True   |
| `script_functions`      |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/scRNA_functions.R                |          | True   |
| `script_preprocess`     |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/seurat_preprocess.Rmd            |          | True   |
| `script_merge`          |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/seurat_merge.Rmd                 |          | True   |
| `script_bc_harmony`     |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/batch_correction_harmony.Rmd     |          | True   |
| `script_bc_rpca`        |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/batch_correction_rpca.Rmd        |          | True   |
| `script_bc_cca`         |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/batch_correction_cca.Rmd         |          | True   |
| `script_liger`          |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/batch_correction_liger.Rmd       |          | True   |
| `script_bc_integration` |             | `string` | /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/SINCLAIR/sinclar-dev-sovacool/bin/batch_correction_integration.Rmd |          | True   |
