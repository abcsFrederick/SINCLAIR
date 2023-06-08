#!/usr/bin/env nextflow
/*
===================================================================
    CCBR/scRNASeq
===================================================================
    GitHub: https://github.com/CCBR/TechDev_scRNASeq_Dev2023
-------------------------------------------------------------------

*/

// Using DSL-2
nextflow.enable.dsl=2

log.info """\
         s c R N A S E Q - N F   P I P E L I N E  
         ===================================
         NF version   : $nextflow.version
         profile      : $workflow.profile
         start time   : $workflow.start
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         outDir       : $params.outdir
         """
         .stripIndent()

//

/*
===================================================================
    Identify workflows
===================================================================
*/
include { PREPROCESS_EXQC                           } from './workflows/pre_process' 
include { GEX_EXQC                                  } from './workflows/gex'
include { ATAC_EXQC                                 } from './workflows/atac'
// include { VDJ_EXQC                            } from '.workflows/sCRNA_vdj'
// include { CITE_EXQC                           } from '.workflows/sCRNA_cite'

//
// WORKFLOW: Run initialization on input samples, check manifests files
//
workflow GEX {
    main:
        PREPROCESS_EXQC ()
        GEX_EXQC (
            PREPROCESS_EXQC.out.ch_meta,
            PREPROCESS_EXQC.out.group_samplesheet,
            PREPROCESS_EXQC.out.h5
        )

}

workflow ATAC {
    main:
        PREPROCESS_EXQC ()
        ATAC_EXQC ()
    emit:
        samplesheet         = PREPROCESS_EXQC.out.samplesheet
}

// //
// // WORKFLOW: Run file sample type 1
// //
// workflow GEX {
//     main:
//         GEX_EXQC()
//     emit:
//         file1                   = GEX_EXQC.out.file1
//         file2                   = GEX_EXQC.out.file2
// }


// //
// // WORKFLOW: Run file sample type 2
// //
// workflow ATAC {
//     main:
//         ATAC_EXQC()
//     emit:
//         file1                   = SUBWK_2_EXQC.out.file1
//         file2                   = SUBWK_2_EXQC.out.file2
// }