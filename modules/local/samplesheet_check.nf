process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    input:
    path (samplesheet)

    output:
    path '*gex_samplesheet.csv'                  , optional:true, emit:gex_samplesheet
    path '*atac_samplesheet.csv'                 , optional:true, emit:atac_samplesheet
    
    script:
    """
    check_samplesheet.py \\
    $samplesheet \\
    project
    """
}