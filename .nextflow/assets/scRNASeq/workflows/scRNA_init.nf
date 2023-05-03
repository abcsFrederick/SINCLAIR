/*
=======================================================================================================
Validate inputs
=======================================================================================================
*/

// Check input path parameters
def chekPathParamList = [ params.input, params.multiqc_config]
for (param in checkPathParamList) { if (param) {file(param, checkIfExists: true)}}

// Check mandatory params
if (params.input) {ch_input = file(params.input)} else {exit 1, 'Input/samplesheet/list not provided!'}

/*
=======================================================================================================
Configs
=======================================================================================================
*/
ch_multiqc_config =                 file("$projectDir/assess/multiqc_config.yaml", checkIfExists:true)

/*
=======================================================================================================
Assign local subworkflows
=======================================================================================================
*/
include { INPUT_CHECK                               } from '../subworkflows/input_check'
include { multiqc_config                            } from '../subworkflows/multiqc'

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SAMPLESHEET_CHECK                         } from '../modules/local/samplesheet_check.nf'
include { MODULE1                                   } from '../modules/local/module1'
include { MODULE2                                   } from '../modules/local/module2'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow SCRNA_INIT_EXQC {
    main:
        // Set output path to relative
        outdir_path = Channel.fromPath(params.outdir,relative:true)

        // Read in samplesheet
        INPUT_CHECK (
            ch_input
        )
}