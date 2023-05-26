process SEURAT_MERGE {
    tag "${gid}"

    input:
    tuple val(gid), val(rdsFiles)
    val(species)
    path(Rlib_dir)
    val(Rpkg_config)    

    output:
    tuple val(gid), path ("*.rds")                  , emit:rds
    
    script:
    def args = task.ext.args ?: ''
    """
    echo $rdsFiles > ${gid}.rds

    Rscript -e 'rmarkdown::render("seurat_merge.Rmd",
        params=list(species=$species,
            rdsFiles=$rdsFiles,
            gid=$gid,
            Rlib_dir=$Rlib_dir,
            Rpkg=$Rpkg,
            testing="N"),
        output_file = "${id}.html"
        output_dir = $task.workDir'
    """

    stub:
    """
    echo $rdsFiles > ${gid}_merged.rds
    """
}