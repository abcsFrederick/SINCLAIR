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
include { INPUT_CHECK_GEX                               } from '../subworkflows/input_check'

/*
=======================================================================================================
Assign Local Modules
=======================================================================================================
*/
include { SAMPLESHEET_CHECK                             } from '../modules/local/samplesheet_check.nf'
include { CELLRANGER_COUNT                              } from '../modules/local/cellranger_count_gex.nf'
include { SEURAT_SINGLE                                 } from '../modules/local/scRNA_single.nf'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow GEX_EXQC {
    main:
        // Set output path to relative, species
        outdir_path = Channel.fromPath(params.outdir,relative:true)
        
        // Read in samplesheet, create metadata
        INPUT_CHECK_GEX (
            ch_input
        )
        
        ch_meta = INPUT_CHECK_GEX.out.gex_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def id = row["sample"]
                    def inDir = row["gex_input_dir"]
                    return [ id, inDir ]
                }
        
        // Run cellranger count
        CELLRANGER_COUNT (
            ch_meta,
            params.genome_dir
        )

        // // Run Seurat for individual samples
        // SEURAT_SINGLE (
        //     ch_meta,
        //     CELLRANGER_COUNT.out.h5,
        //     params.species,
        //     params.group_tab,
        //     params.Rlib_dir
        // )

    emit:
        samplesheet        = INPUT_CHECK_GEX.out.gex_samplesheet
}