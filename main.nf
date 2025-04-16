#!/usr/bin/env nextflow
/*
===================================================================
    CCBR/SINCLAIR
===================================================================
    GitHub: https://github.com/CCBR/SINCLAIR
-------------------------------------------------------------------
*/

// Using DSL-2
nextflow.enable.dsl=2

/*
===================================================================
    Import workflows
===================================================================
*/
include { PREPROCESS_EXQC                           } from './workflows/pre_process'
include { GEX_EXQC                                  } from './workflows/gex'
include { ATAC_EXQC                                 } from './workflows/atac'
// include { VDJ_EXQC                            } from '.workflows/sCRNA_vdj'
// include { CITE_EXQC                           } from '.workflows/sCRNA_cite'

// Plugins
include { validateParameters } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'


/*
===================================================================
    Set handlers
===================================================================
*/
workflow.onComplete {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}

/*
===================================================================
    Define workflows
===================================================================
*/
//
// WORKFLOW: Run initialization on input samples, check manifests files
//
workflow {
    main:
        log.info """\
                SINCLAIR $workflow.manifest.version
                ===================================
                NF version   : $nextflow.version
                runName      : $workflow.runName
                username     : $workflow.userName
                configs      : $workflow.configFiles
                profile      : $workflow.profile
                cmd line     : $workflow.commandLine
                start time   : $workflow.start
                projectDir   : $workflow.projectDir
                launchDir    : $workflow.launchDir
                workDir      : $workflow.workDir
                homeDir      : $workflow.homeDir
                outDir       : $params.outdir
                """
                .stripIndent()

        validateParameters()
        PREPROCESS_EXQC ()
        GEX_EXQC (
            PREPROCESS_EXQC.out.ch_fqdir_h5,
            PREPROCESS_EXQC.out.group_samplesheet,
        )

}
