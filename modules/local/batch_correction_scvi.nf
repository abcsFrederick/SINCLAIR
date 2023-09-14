process PLACEHOLDER_SCVI {
    tag "${gid}"
    conda 'scvi-env'

    input:
    tuple val(gid), path(rdsFiles)
    val(species)
    val(Rlib_dir)
    path(Rpkg_config)
    path(rmd)
    path(scRNA_functions)  

    output:
    tuple val(gid), path ("*_seurat_merged.rds")                 , emit:rds
    tuple val(gid), path ("*_seurat_merged.pdf")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    # need to source the pipeliner conda env
    . "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
    
    Rscript -e 'rmarkdown::render("${rmd}",
        params=list(species="$species",
            rdsFiles="$rdsFiles",
            gid="$gid",
            Rlib_dir="$Rlib_dir",
            Rpkg_config="$Rpkg_config",
            scRNA_functions="$scRNA_functions",
            testing="N"),
        output_file = "${gid}_seurat_merged.pdf")'
    """

    stub:
    """
    echo $rdsFiles > ${gid}_merged.rds
    touch ${gid}.pdf
    """
}