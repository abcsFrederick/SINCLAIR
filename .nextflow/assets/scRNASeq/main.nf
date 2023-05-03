#!/usr/bin/env nextflow
/*
===================================================================
    CCBR/scRNASeq
===================================================================
    GitHub: https://github.com/CCBR/TechDev_scRNASeq_Dev2023
-------------------------------------------------------------------

*/

/*
===================================================================
    Identify workflows
===================================================================
*/
include { SCRNA_INIT_EXQC                } from '.workflows/scRNA_init'
include { SUBWK_1_EXQC                   } from '.workflows/scRNA_sub1'
include { SUBWK_2_EXQC                   } from '.workflows/sCRNA_sub2'

//
// WORKFLOW: Run initialization on input samples, check manifests files
//
workflow SCRNA_INIT_EXQC {
    main:
        SCRNA_INIT_EXQC()
    emit:
        file1                   = SCRNA_INIT_EXQC.out.file1
        file2                   = SCRNA_INIT_EXQC.out.file2
}

//
// WORKFLOW: Run file sample type 1
//
workflow SUBWK_1 {
    main:
        SUBWK_1_EXQC()
    emit:
        file1                   = SUBWK_1_EXQC.out.file1
        file2                   = SUBWK_1_EXQC.out.file2
}


//
// WORKFLOW: Run file sample type 2
//
workflow SUBWK_2 {
    main:
        SUBWK_2_EXQC()
    emit:
        file1                   = SUBWK_2_EXQC.out.file1
        file2                   = SUBWK_2_EXQC.out.file2
}