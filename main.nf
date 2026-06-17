#!/usr/bin/env nextflow
/*
 * ========================================================================================
 *  pals_celloracle — Nextflow (DSL2) port of the CellOracle scvelo-pancreas tutorial
 * ========================================================================================
 *  Github   : https://github.com/BCI-KRP/pals_celloracle
 *  Tutorial : https://github.com/morris-lab/CellOracle
 *  Nextflow : https://www.nextflow.io/docs/latest/index.html
 * ----------------------------------------------------------------------------------------
 *  Usage
 *    nextflow run main.nf -profile local                 # laptop/workstation
 *    nextflow run main.nf -profile apocrita              # QMUL Apocrita SLURM cluster
 *    nextflow run main.nf --help                         # show this help text
 * ----------------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------------------
//  Shared CLI helpers (help banner, run summary)
// ----------------------------------------------------------------------------------------
include { helpMessage; logRunSummary; containerPullHint } from './lib/help.nf'

// ----------------------------------------------------------------------------------------
//  Pipeline logic — one sub-workflow, one line.
// ----------------------------------------------------------------------------------------
include { CELLORACLE } from './subworkflows/local/celloracle.nf'

// ----------------------------------------------------------------------------------------
//  Main workflow
// ----------------------------------------------------------------------------------------
workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    logRunSummary()
    containerPullHint()

    CELLORACLE()
}
