#!/bin/bash

module load cellranger

cellranger count --id <outputDirName> \
        --fastqs <path/to/fastq/dir> \
        --transcriptome=<path/to/genome> \
        --localcores=$SLURM_CPUS_PER_TASK \
        --localmem=34 \
        --jobmode=slurm --maxjobs=20

#This command is run for each sample
#FASTQ directory should contain paired-end read files with some unknown reads in the format R1, R2, or I1. There are normally 3 or 4 sets per sample
#The genome reference directory should contain 4 directories (fasta, genes, pickle, star) and a reference.json. Refer to /fdb/cellranger/refdata-cellranger-3.0.0/mm10 for an example
