process SEURAT_PREPROCESS {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    path(h5)
    val(species)
    path(Rlib_dir)
    path(Rpkg)

    output:
    tuple val(id), path ("*.rds")                  , emit:rds
    tuple val(id), path ("*.html")                 , emit:logs
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("seurat_preprocess.Rmd",
        params=list(species=$species,
            id=$id,
            h5=$h5,
            Rlib_dir=$Rlib_dir,
            Rpkg=$Rpkg,
            testing="N"),
        output_file = "${id}.html"
        output_dir = $task.workDir'
    """

    stub:
    """
    touch ${id}_seurat_object.rds
    touch ${id}.html
    """
}