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
include { SCRNA_INIT_EXQC                     } from '.workflows/scRNA_init'
include { GEX_EXQC                            } from '.workflows/scRNA_gex'
include { ATAC_EXQC                           } from '.workflows/sCRNA_atac'
include { VDJ_EXQC                            } from '.workflows/sCRNA_vdj'
include { CITE_EXQC                           } from '.workflows/sCRNA_cite'

//
// WORKFLOW: Run initialization on input samples, check manifests files
//
workflow SCRNA_INIT {
    main:
        SCRNA_INIT_EXQC()
    emit:
        file1                   = SCRNA_INIT_EXQC.out.file1
        file2                   = SCRNA_INIT_EXQC.out.file2
}

//
// WORKFLOW: Run file sample type 1
//
workflow GEX {
    main:
        GEX_EXQC()
    emit:
        file1                   = SUBWK_1_EXQC.out.file1
        file2                   = SUBWK_1_EXQC.out.file2
}


//
// WORKFLOW: Run file sample type 2
//
workflow ATAC {
    main:
        SUBWK_2_EXQC()
    emit:
        file1                   = SUBWK_2_EXQC.out.file1
        file2                   = SUBWK_2_EXQC.out.file2
}