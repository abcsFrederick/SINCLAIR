process SEURAT_MERGE {
    tag "${gid}"

    input:
    tuple val(gid), val(rdsFiles)
    val(species)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)  

    output:
    tuple val(gid), path ("*.rds")                , emit:rds
    tuple val(id), path ("*.pdf")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(species=$species,
            rdsFiles=$rdsFiles,
            gid=$gid,
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}.pdf")'
    """

    stub:
    """
    echo $rdsFiles > ${gid}_merged.rds
    touch ${gid}.pdf
    """
}