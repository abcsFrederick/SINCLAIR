process SEURAT_SINGLE {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    path(h5)
    val(species)
    path(Rlib_dir)

    output:
    path '*.rds'                  , emit:rds
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript scRNA.R \
        $species \
        $id \
        $h5 \
        $Rlib_dir
    """

    stub:
    """
    touch ${id}_seurat_object.rds
    """
}