# Expected Outputs

The following folders and files are created under the working directory SINCLAIR was initialized in:

- **batch_correct**: contains the various RDS files for batch correction methods (CCA, HARMONY, scVI, RPCA, LIGER) and batch correction report (HTML)

A detailed walkthrough for each method can be found in:

[CCA](https://direct.mit.edu/neco/article/16/12/2639/6880/Canonical-Correlation-Analysis-An-Overview-with)

[HARMONY](https://portals.broadinstitute.org/harmony/advanced.html)

[scVCI](https://docs.scvi-tools.org/en/0.6.8/tutorials/basic_tutorial.html)

[RPCA](https://satijalab.org/seurat/articles/integration_rpca)

[LIGER](https://welch-lab.github.io/liger/index.html)

- **cellranger_counts**: contains the h5 files, if `cellranger count` is deployed
- **pipeline_info**: contains execution_reports, execution_trace and pipeline_dag files from NextFlow
- **samplesheets**: contains the manifests used to identify samples, contrasts, and sample:contrast groupings
- **seurat**: contains multiple seurat-generated directories:
  - preprocess: contains sample level data; both RDS and PDF files from pre-processing
  - merge: contrast grouped, sample level data; both RDS and PDf files of pre-processed, merging

```
results
├── batch_correct
│   ├── group1-group2_batch_correction_cca.rds
│   ├── group1-group2_batch_correction_harmony.rds
│   ├── group1-group2_batch_correction_integration.html
│   ├── group1-group2_batch_correction_liger.rds
│   └── group1-group2_batch_correction_rpca.rds
├── cellranger_counts
│   ├── sample1
│   │   └── outs
│   │       └── filtered_feature_bc_matrix.h5
│   ├── ...
│   └── sampleN
│       └── outs
│           └── filtered_feature_bc_matrix.h5
├── samplesheets
│   ├── project_contrast_samplesheet.csv
│   ├── project_gex_samplesheet.csv
│   └── project_groups_samplesheet.csv
└── seurat
    ├── merge
    │   └── group1-group2_seurat_merged.rds
    └── preprocess
        ├── sample1_seurat_preprocess.rds
        ├── ...
        └── sampleN_seurat_preprocess.rds
```

The following table contains a description of each file and their potential use(s) cases for single-cell analysis.

| File/Folder Path                                                | Description                                                                                                                                                                                                                                             | Use Case(s)                                                                                                                                                                                      |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `batch_correct/group1-group2_batch_correction_<method>.rds`     | Seurat object after batch correction using one of the available methods.                                                                                                                                                                                | For downstream analyses (e.g. Differential Gene Expression, Pathway) that requires batch corrected data. Consult the batch_correction_integration.html for which batch corrected dataset to use. |
| `batch_correct/group1-group2_batch_correction_integration.html` | HTML report summarizing integration diagnostics and visualizations. Includes dimension reduction plots (e.g., UMAP) of the integrated data annotated by cell types, and silhouette plots to assess the quality of integration across different methods. | Assess which batch correction method(s) best resolve technical differences across samples and can be used for downstream analyses                                                                |
| `cellranger_counts/sampleN/outs/filtered_feature_bc_matrix.h5`  | Filtered gene expression matrix from **Cell Ranger**; used for downstream analysis.                                                                                                                                                                     |                                                                                                                                                                                                  |
| `samplesheets/project_gex_samplesheet.csv`                      | Metadata mapping samples for gene expression input (e.g., sample names, paths).                                                                                                                                                                         |                                                                                                                                                                                                  |
| `samplesheets/project_groups_samplesheet.csv`                   | Defines sample-to-group mapping for batch correction or merging.                                                                                                                                                                                        |                                                                                                                                                                                                  |
| `samplesheets/project_contrast_samplesheet.csv`                 | Specifies contrast groups for differential expression testing.                                                                                                                                                                                          |                                                                                                                                                                                                  |
| `seurat/merge/group1-group2_seurat_merged.rds`                  | Merged Seurat object combining multiple samples or groups.                                                                                                                                                                                              |                                                                                                                                                                                                  |
| `seurat/preprocess/sampleN_seurat_preprocess.rds`               | Seurat object after QC and normalization for an individual sample.                                                                                                                                                                                      |                                                                                                                                                                                                  |
