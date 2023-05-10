/*
=======================================================================================================
Validate inputs
=======================================================================================================
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//input on command line
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet/list not specified!' }

/*
=======================================================================================================
Assign local subworkflows
=======================================================================================================
*/
include { INPUT_CHECK                               } from '../subworkflows/input_check'

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SAMPLESHEET_CHECK                         } from '../modules/local/samplesheet_check.nf'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow GEX_EXQC {
    main:
        // Set output path to relative
        outdir_path = Channel.fromPath(params.outdir,relative:true)

        // Read in samplesheet
        INPUT_CHECK (
            ch_input
        )
    emit:
        samplesheet        = INPUT_CHECK.out.samplesheet
}