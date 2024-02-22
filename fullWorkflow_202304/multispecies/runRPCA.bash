#!/bin/bash

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
#conda init scRNA4
conda activate scRNA4

module load R/4.2
#module load python

Rscript workflow/scripts/rpca.R filtered integration/allSamples_rpca.rds both all 0.2,0.4,0.6,0.8,1.0,1.2,1.4
