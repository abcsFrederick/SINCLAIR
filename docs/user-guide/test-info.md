# Setting up a test run

A test run of the pipeline can be initiated through the following command:

```
sinclair run --mode slurm -profile test
```

# About the data

**_Overview_**

**Species:** Homo sapiens

**Genome Reference:** Human genome GRCh38

**Data Type:** Single-cell gene expression (GEX)

**Sequencing Platform:** 10x Genomics Chromium, 3' v3.1 chemistry

**Analysis Software** Cell Ranger v6.1.0

## **INFO:**

**Donor Information:** healthy female.

**Isolation protocol:** CG000392 RevA: Isolation of Leukocytes, Bone Marrow and Peripheral Blood Mononuclear Cells for Single Cell RNA Sequencing - Whole Blood Lysis for Granuloctyes track.

Whole transcriptome/Gene Expression libraries were generated as described in the Chromium Next GEM Single Cell 3' Reagent Kits v3.1 (Dual Index) User Guide (CG000204 Rev D).

**TUTORIAL:** https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/tutorials/neutrophils

**SOURCE:** https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard

# Sample Structure and Setup

The original dataset was downloaded from 10x Genomics in FASTQ format and unpacked from WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs.tar. It is divided into four sets of FASTQ files:

## WB_1

- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_I1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_I2_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_R1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L001_R2_001.fastq.gz

## WB_2

- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_I1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_I2_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_R1_001.fastq.gz
- WB_Lysis_Granulocytes_3p_Introns_8kCells_S1_L002_R2_001.fastq.gz

## WB_3 and WB_4 - Simulated Replicates

To simulate technical and biological replicates, and to assess the effects of clustering, batch correction, and downstream steps (e.g. Differential Gene Expression), WB_3 and WB_4 samples were copied from WB_1 and WB_2 respectively.

- WB_1 --> WB_3
- WB_2 --> WB_4

## Subsampling Details

Each group was further downsampled to make the dataset lightweight and suitable for rapid testing.

- WB_1 --> sample1,sample2 --> s100
- WB_2 --> sample3,sample4 --> s101
- WB_3 --> sample5,sample6 --> s102
- WB_4 --> sample7,sample8 --> s103

**_Test Grouping Setup_**

For testing group-based logic

# Single Sample

<<<<<<< Updated upstream

# for single sample

- group1

# for multiple samples

=======

- group1

# Multiple Samples

> > > > > > > Stashed changes

- group1,group2
- group1,group2,group3

## Expected Output

The SINCLAIR test run will produce:

Filtered gene-barcode matrix

Batch corrected UMAP/t-SNE dimensionality reduction plots, annotated by celltypes

- Technical replicates (WB_1 + WB_3 and WB_2 + WB_4) should colocalize in the same regions of the dimension reduction plots.

Cluster assignments, Silhouette scores quantifying the degree of integration occurring in each cluster, per batch correction method.

- Clusters should be separated on the basis of celltype annotations rather than batch origins.

## Issues and Troubleshooting

Come across a **bug**? Open an [issue](https://github.com/CCBR/SINCLAIR/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/SINCLAIR/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](../contributing.md).

**General Inquiries and Collaboration:** Please contact the CCBR Pipeliner team at [CCBR_Pipeliner@mail.nih.gov](mailto:CCBR_Pipeliner@mail.nih.gov).
