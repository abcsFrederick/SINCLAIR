process BATCHC_INT {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    path(singleRDS)
    val(species)
    val(int_params)
    path(contrast)
    path(Rlib_dir)

    output:
    path '*.rds'                  , emit:rds
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript integrateBatches.R \
		$singleRDS \
        seuratIntegrated.rds  \
        $species \
        $contrasts \
        $int_params \
        $Rlib_dir \
    """
}