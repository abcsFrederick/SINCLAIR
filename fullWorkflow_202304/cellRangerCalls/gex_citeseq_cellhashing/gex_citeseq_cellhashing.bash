#!/bin/bash

module load cellranger

cellranger count --id=<outputDirName> \
--transcriptome=<path/to/genome/> \
--libraries=<sample_library_manifest.csv> \
--feature-ref=features.csv \
--localcores=$SLURM_CPUS_PER_TASK \
--localmem=34

#Primary difference between this and basic GEX is the --libraries and --feature-ref flags. The sample library manifest should be sample specific, but the features.csv file is more likely to be generic
