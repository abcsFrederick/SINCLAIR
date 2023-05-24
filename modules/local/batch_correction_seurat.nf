process SEURAT_INTEGRATION {
    tag "${contrast_id}"

    input:
    tuple val(contrast_id), val(ct1), val(ct2)
    val(int_params)
    val(species)
    path(Rlib_dir)
    path(Rpkg)

    output:
    path '*.rds'                  , emit:rds
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript integrateBatches.R \
        $species \
        $ct1,$ct2 \
        $int_params\
        $Rlib_dir \
        $Rpkg
    """

    stub:
    """
    touch ${contrast_id}_seurat_object.rds
    """
}