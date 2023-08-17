# Getting started:

The 3 following `.csv` files need to be modified to reflect the data being analyzed.

### `input_manifest.csv`
This file is used when starting with the CellRanger alignment. It resembles the library.csv files required for initiating a [10X CellRanger multi run](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#libraries). It is a comma-separated file with the following fields:

* **masterID**: A parent identifier that can be used to group technical replicates of the same sample.
* **uniqueID**: A distinct identifier used for each individual replicate.
* **groupID**: Identifier used for experimental groups (e.g. Treatment and Control).
* **dataType**: Identifier for single cell data type. Currently supports only `gex`. Upcoming options will include `citeseq`, `vdj`, and `atac`.
* **inputDir**: The full path to the `fastq.gz` files used for CellRanger alignment.

### `input_manifest_cellranger.csv`
This file is used when starting with aligned files, specifically `.h5` compressed files that are produced by 10X CellRanger. It has the same fields as the `input_manifest.csv`, with the following modification:

* **inputDir**: The full path to the `.h5` files produced by CellRanger alignment. This typically takes the form of `/<path>/<to>/<sample_output>/outs`.

### `contrast_manifest.csv`
This file is used for determining the contrasts to be compared, when groups are defined by samples. These contrasts are then used in the sample combination steps of `merge` and `batch_correct`. Multiple contrasts can be submitted, and samples not within the contrast groups will be excluded from that particular contrast. It is a comma-separated file in the following format, where each line represents a contrast:

`contrast1,contrast2,contrast3,...,contrast_n`
`group1,group2`
`group1,group2,group3`
`group1,group3`

In this case, 3 contrasts will be generated:

* group1-group2
* group1-group2-group3
* group1-group3
