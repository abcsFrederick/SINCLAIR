# 4. Expected Outputs
The following directories are created under the WORKDIR directory:

- cellranger_counts: contains the h5 files, if `cellranger count` is deployed
- pipeline_info: contains execution_reports, execution_trace and pipeline_dag files from NextFlow
- samplesheets: contains the manifests used to identify samples, contrasts, and sample:contrast groupings
- seurat: contains multiple seurat-generated directories:
    - preprocess: contains sample level data; both RDS and PDF files from pre-processing
    - merge: contrast grouped, sample level data; both RDS and PDf files of pre-processed, merging
    - batch_correction: contrast grouped level data; contains directories for each batch_correction
        - harmony
        - cca
        - rpca

```
cellranger_counts
│   ├── sample1
│   │   └── outs
seurat
    ├── merge
    │   ├── group1-group2-group3_seurat_merged.pdf
    │   ├── group1-group2-group3_seurat_merged.rds
    │   ├── group1-group2_seurat_merged.pdf
    │   └── group1-group2_seurat_merged.rds
    └── preprocess
        ├── sample1_seurat_preprocess.pdf
        ├── sample1_seurat_preprocess.rds
        ├── sample2_seurat_preprocess.pdf
```