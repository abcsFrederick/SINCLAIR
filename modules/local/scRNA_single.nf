process SEURAT_SINGLE {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    path(h5)
    val(species)
    path(group_tab)
    path(Rlib_dir)

    output:
    path '*.rds'                  , emit:rds
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript scRNA.R \
        $h5 \
        $species \
        ${id}.rds \
        $group_tab \
        $Rlib_dir
    """
}