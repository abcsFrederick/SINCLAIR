process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    input:
    path samplesheet

    outptu:
    path '*.valid.csv', emit:csv
    
    script:
    """
    check_samplesheet.py \\
    $samplesheet \\
    samplesheet.valid.csv
    """
}