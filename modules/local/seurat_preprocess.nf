process SEURAT_PREPROCESS {
    tag "${id}"
    label 'process_medium'

    container "${params.container_seurat}"

    input:
    tuple val(id), val(inDir), path(h5), path(celldex_path)
    val(species)
    val(qc_filtering)
    val(nCount_RNA_max)
    val(nCount_RNA_min)
    val(nFeature_RNA_max)
    val(nFeature_RNA_min)
    val(percent_mt_max)
    val(percent_mt_min)
    val(run_doublet_finder)
    val(npcs)
    path(rmd)
    path(scRNA_functions)

    output:
    tuple val(id), path ("*.rds")                 , emit:rds
    tuple val(id), path ("*.html")                 , emit:logs

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript --vanilla
    rmarkdown::render("${rmd}",
        params=list(
            species="$species",
            sampleid="$id",
            h5="$h5",
            qc_filtering="$qc_filtering",
            nCount_RNA_max=$nCount_RNA_max,
            nCount_RNA_min=$nCount_RNA_min,
            nFeature_RNA_max=$nFeature_RNA_max,
            nFeature_RNA_min=$nFeature_RNA_min,
            percent_mt_max=$percent_mt_max,
            percent_mt_min=$percent_mt_min,
            run_doublet_finder="$run_doublet_finder",
            npcs=$npcs,
            scRNA_functions="$scRNA_functions",
            celldex_cache="$celldex_path"
            ),
        output_file = "${id}_seurat_preprocess.html"
    )
    """

    stub:
    """
    touch ${id}_seurat_preprocess.rds
    touch ${id}_seurat_preprocess.html
    """
}
