process SEURAT_PREPROCESS {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    path(h5)
    val(species)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(id), path ("*.rds")                 , emit:rds
    tuple val(id), path ("*.pdf")                 , emit:logs
    
    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(species="$species",
            sampleid="$id",
            h5="$h5",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${id}_seurat_preprocess.pdf")'
    """

    stub:
    """
    touch ${id}_seurat_preprocess.rds
    touch ${id}_seurat_preprocess.pdf
    """
}