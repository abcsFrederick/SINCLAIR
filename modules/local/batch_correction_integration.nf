process BATCH_CORRECT_INTEGRATION {
    tag "${gid}"
    label 'process_high'

    input:
    tuple val(gid), path(mergedObj)
    tuple val(gid), val(rds_h)
    tuple val(gid), val(rds_r)
    tuple val(gid), val(rds_c)
    tuple val(gid), val(rds_s)
    val(species)
    val(npcs)
    val(resolution_list)
    path(Rlib_dir)
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
            mergedObj="$rds_m",
            harmonyObj="$rds_h",
            ccaObj="$rds_c",
            rpcaObj="$rds_r",
            scviObj="$rds_s",
            npcs="$npcs",
            resolution_list="$resolution_list",
            citseq="",
            annot="",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_batch_correction_integration.pdf")'
    """

    stub:
    """
    touch ${gid}_batch_correction_integration.rds
    touch ${gid}_batch_correction_integration.pdf
    """
}
