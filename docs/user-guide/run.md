# 3. Running the Pipeline

## 3.1 Pipeline Overview
The Nextflow workflow can be run as follows:
```
nextflow run main.nf \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \
    $args
```

## 3.2 Commands explained
The following explains each of the command options:

- entry: accepts the datatype to be used; IE gex
- profile: how to run the processes; IE biowulf singularity, docker
- input: input_manifest.csv location
- contrast: contrast_manifest.csv location
- outdir: complete path to the output dir
- species: species to be used
- args: any additional arguments; IE --stub-run

## 3.3 Typical Workflow
A typical command workflow, running the pipeline for the first time locally, is as follows:
```
nextflow run main.nf \
    -entry gex \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --species hg19
```

A typical command workflow, running the pipeline for a repeated time locally, is as follows:
```
nextflow run main.nf -resume \
    -entry gex \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --species hg19
```

A typical command workflow, running the pipeline in a `dryrun mode`, is as follows:
```
nextflow run main.nf \
    -entry gex \
    -profile biowulf \
    --run_cellranger Y \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /path/to/scRNA_test \
    --species hg19 \
    --stub-run
```

Alternatively a script was created to run the pipeline, which takes the following flags:
- species: hg19
- datatype: GEX
- resume: Y or N
- outDir: /path/to/output/dir
- args: --stub-run, --cellranger Y

Examples:
```
# run GEX on test data, for the first time
sh run_scRNA.sh hg19 GEX N /path/to/output/dir

# run GEX on test data, as a re-run
sh run_scRNA.sh hg19 GEX Y /path/to/output/dir

# run GEX on test data, as a dry-run
sh run_scRNA.sh hg19 GEX N /path/to/output/dir -stub-run 
```