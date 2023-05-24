process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    input:
    path (samplesheet)
    path (contrast_samplesheet)

    output:
    path '*gex_samplesheet.csv'                      , optional:true, emit:gex_samplesheet
    path '*atac_samplesheet.csv'                     , optional:true, emit:atac_samplesheet
    path '*contrast_samplesheet.csv'                 , emit:contrast_samplesheet
    path '*group_samplesheet.csv'                    , emit:group_samplesheet
    
    script:
    """
    check_samplesheet.py \\
    $samplesheet \\
    $contrast_samplesheet \\
    project
    """

    stub:
    """
    check_samplesheet.py \\
    $samplesheet \\
    $contrast_samplesheet \\
    project
    """
}