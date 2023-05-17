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
if (params.input)    { ch_input    = file(params.input)    } else { exit 1, 'Input samplesheet not specified!' }
if (params.contrast) { ch_contrast = file(params.contrast) } else { exit 1, 'Contrast samplesheet not specified!' }

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
        
        // Read in samplesheet, contrast manifest
        INPUT_CHECK_GEX (
            ch_input,
            ch_contrast
        )
        
        // create metadata
        ch_meta = INPUT_CHECK_GEX.out.gex_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def id = row["sample"]
                    def inDir = row["gex_input_dir"]
                    return [ id, inDir ]
                }

        ch_contrast = INPUT_CHECK_GEX.out.contrast_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def contrast1 = row["contrast1"]
                    def contrast2 = row["contrast2"]
                    return [ contrast1, contrast2 ]
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
        //     params.Rlib_dir
        // )

        // // Run Merge TBD from AA
        // SEURAT_MERGE (
        //     ch_meta,
        //     CELLRANGER_COUNT.out.h5,
        //     params.species,
        //     params.Rlib_dir
        // )

        // // Run batch corrections
        // BATCHC_INT (
        //     ch_contrast,
        //     SEURAT_MERGE.out.rds,
        //     params.species,
        //     params.seurat_integration,
        //     params.Rlib_dir
        // )
        
        // // Run batch corrections
        // BATCHC_RPCA (
        //     ch_meta,
        //     CELLRANGER_COUNT.out.h5,
        //     params.species,
        //     params.Rlib_dir
        // )

        // // Run batch corrections
        // BATCHC_HARMONY (
        //     ch_meta,
        //     CELLRANGER_COUNT.out.h5,
        //     params.species,
        //     params.Rlib_dir
        // )
    emit:
        samplesheet        = INPUT_CHECK_GEX.out.gex_samplesheet
}