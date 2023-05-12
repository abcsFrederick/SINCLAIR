process CELLRANGER_COUNT {
    tag "${id}"

    input:
    tuple val(id), val(inDir)
    val(genome_dir)

    output:
    path '*/outs'                  , emit:cr_outDir
    
    script:
    def args = task.ext.args ?: ''
    def localmem = task.memory.toGiga()
    """
    cellranger count \
        --id $id \
        --fastqs $inDir \
        --transcriptome=${genome_dir} \
        --localcores=$task.cpus \
        --localmem=$localmem \
        $args
    """
}