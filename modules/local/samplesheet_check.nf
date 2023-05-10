process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    input:
    path (samplesheet)

    output:
    path '*_samplesheet.csv', emit:csv
    
    script:
    """
    check_samplesheet.py \\
    $samplesheet \\
    project
    """
}