# SINCLAIR - **SIN**gle **C**el**L** **A**nalys**I**s **R**esource

## Table of contents
- [SINCLAIR - **SIN**gle **C**el**L** **A**nalys**I**s **R**esource](#sinclair---single-cell-analysis-resource)
  - [Table of Contents](#table-of-contents)
  - [1. Introduction](#1-Introduction)
  - [2. Overview](#2-Overview)
  - [3. Documentation](#3-Documentation)
  - [4. Feedback](#4-Feedback)

### 1. Introduction
The [SINCLAIR - **SIN**gle **C**el**L** **A**nalys**I**s **R**esource](#sinclair---single-cell-analysis-resource) was developed by the CCR Collaborative Bioinformatics Resource as an open-source, reproducible solution for multiple single cell next-generation modalities. It has been developed soley on [Biowulf](https://hpc.nih.gov/) using [Nextflow](https://www.nextflow.io/).

### 2. Overview
The pipeline currently begins with either sample FASTQ file or h5 Aligned reads, completing per sample quality control, and per-contrast integration. Quality control reports are generated, as are per-contrast integration reports.

![Single cell RNA-Seq GEX pipeline](./resources/scRNA.svg) <sup>**Overview of Single Cell RNASeq Gene Expression Pipeline**</sup>

### 3. Documentation
Please view the repositories [documentation](https://symmetrical-adventure-ovjq9gl.pages.github.io/) for full details on deploying the pipeline, features, testing, and expected outputs.

### 4. Feedback
For comments/suggestions/advice please reach out to [Samantha Chill](mailto:samantha.sevilla@nih.gov) or [Nathan Wong](mailto:nathan.wong@nih.gov).
