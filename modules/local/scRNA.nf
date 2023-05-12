rule perSample:
	input:
		join(workpath,"cellrangerOut","{name}","outs","filtered_feature_bc_matrix.h5")
	output:
		join(workpath,"seuratOut","{name}.rds")
	params:
		specie=specie
	shell: """
. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
#conda init scRNA4
conda activate scRNA4
module load R/4.2.0;
Rscript workflow/scripts/scRNA.R {input} {params.specie}  {output}
