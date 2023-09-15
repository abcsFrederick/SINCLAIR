process BATCH_CORRECT_SCVI {
    tag "${gid}"
    label 'process_high'

    input:
    tuple val(gid), path(mergedObj)
    val(species)
    val(npcs)
    val(resolution_list)
    val(conda_path)
    val(python_path)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)  

    output:
    tuple val(gid), path ("*.rds")                 , emit:rds
    tuple val(gid), path ("*.pdf")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(gid="$gid",
            mergedObj="$mergedObj",
            species="$species",
            npcs="$npcs",
            resolution_list="$resolution_list",
            conda_path="$conda_path",
            python_path="$python_path",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_batch_correction_scvi.pdf")'
    """

    stub:
    """
    touch ${gid}_batch_correction_scvi.rds
    touch ${gid}_batch_correction_scvi.pdf
    """
}