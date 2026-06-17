/*
 * =====================================================================================
 * subworkflows/local/celloracle.nf
 * -------------------------------------------------------------------------------------
 * Chains the five atomic processes together and exposes a single
 * reusable workflow entry point.
 * =====================================================================================
 */

include { CO_DATA_PREP     } from '../../modules/celloracle/data_prep/main.nf'
include { CO_PSEUDOTIME    } from '../../modules/celloracle/pseudotime/main.nf'
include { CO_GRN_INFERENCE } from '../../modules/celloracle/grn_inference/main.nf'
include { CO_PERTURBATION  } from '../../modules/celloracle/perturbation/main.nf'
include { CO_GRADIENT_DEV  } from '../../modules/celloracle/gradient_dev/main.nf'

workflow CELLORACLE {

    // -------- Step 1: preprocessing --------
    if (params.run_preprocess) {
        CO_DATA_PREP(params.dataset)
        ch_adata_prep = CO_DATA_PREP.out.adata
    } else {
        def cached = file("${params.outdir}/01_preprocess/adata_preprocessed.pkl")
        if (!cached.exists()) {
            error "Preprocess was skipped (--run_preprocess false) but cached file is missing: ${cached}"
        }
        ch_adata_prep = channel.fromPath(cached)
    }

    // -------- Step 2: pseudotime --------
    if (params.run_pseudotime) {
        CO_PSEUDOTIME(ch_adata_prep)
        ch_adata_pt = CO_PSEUDOTIME.out.adata
    } else {
        ch_adata_pt = ch_adata_prep
    }

    // -------- Step 3: GRN inference --------
    if (params.run_grn) {
        CO_GRN_INFERENCE(ch_adata_pt)
        ch_oracle = CO_GRN_INFERENCE.out.oracle
        ch_links  = CO_GRN_INFERENCE.out.links
    } else {
        def cached_oracle = file("${params.outdir}/03_grn/oracle_object.pkl")
        def cached_links  = file("${params.outdir}/03_grn/links_object.pkl")
        
        if (!cached_oracle.exists() || !cached_links.exists()) {
            error "GRN was skipped (--run_grn false) but cached files are missing in ${params.outdir}/03_grn"
        }
        ch_oracle = channel.fromPath(cached_oracle)
        ch_links  = channel.fromPath(cached_links)
    }

    // -------- Step 4: in-silico perturbation --------
    if (params.run_perturbation) {
        CO_PERTURBATION(ch_adata_pt, ch_oracle, ch_links)
        ch_oracle_pert = CO_PERTURBATION.out.oracle
    } else {
        ch_oracle_pert = ch_oracle
    }

    // -------- Step 5: gradient + dev module --------
    if (params.run_gradient) {
        CO_GRADIENT_DEV(ch_adata_pt, ch_oracle_pert)
    }
}