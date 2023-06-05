process BATCH_CORRECT_HARMONY {
    tag "${gid}"

    input:
    tuple val(gid), path(mergedObj)
    val(species)
    val(ncps)
    val(resolution_list)
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
        params=list(gid="$gid",
            mergedObj="$mergedObj",
            species="$species",
            ncps="$ncps",
            resolution_list="$resolution_list",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${id}_batch_correction_harmony.pdf")'
    """

    stub:
    """
    touch ${id}_batch_correction_harmony.rds
    touch ${id}_batch_correction_harmony.pdf
    """
}