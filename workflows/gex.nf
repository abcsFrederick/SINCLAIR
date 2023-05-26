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
include { SEURAT_MERGE                                 } from '../modules/local/scRNA_merge.nf'
// include { BATCH_CORRECT_INT                             } from '../modules/local/batch_correction_seurat.nf'

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
        // creates GEX mapped input of sample:gex_input_dir
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

        // Run Seurat for individual samples
        SEURAT_SINGLE (
            ch_meta,
            CELLRANGER_COUNT.out.h5,
            params.species,
            params.Rlib_dir,
            params.Rpkg
        )

        // creates metadata
        // creates tuple of group_id_key: all RDS files from INPUT_CHECK_GEX
        //// which matches those samples
        //// example: INPUT_CHECK_GEX.out.group_samplesheet; SEURAT_SINGLE.out.rds
        //// example: sample1:group1_group2; sample1:sample1.rds
        //// example: sample1:group1_group2_group3; sample1:sample1.rds
        //// example: sample2:group1_group2; sample2:sample2.rds
        //// output: group1_group2: [sample1.rds,sample2.rds]
        //// output: group1_group2_group3: [sample1.rds]
        ch_groups = INPUT_CHECK_GEX.out.group_samplesheet
            .splitCsv( header:true, sep:',', strip:true )
            .map { row ->
                    def key = row["keyid"]
                    def sample = row["sampleid"]
                    return [sample, key]
                }
            .combine(SEURAT_SINGLE.out.rds, by: 0)
            .map { sample, key, rds_file -> tuple( key, rds_file ) }
            .groupTuple()
            .view()


        SEURAT_MERGE (
            ch_groups,
            // INPUT_CHECK_GEX.out.contrast_samplesheet,
            // params.species,
            // params.Rlib_dir,
            // params.Rpkg
        )

        // Run batch corrections
        // BATCH_CORRECT_INT (
        //     ch_contrast, SEURAT_SINGLE.out.rds.collect(),
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