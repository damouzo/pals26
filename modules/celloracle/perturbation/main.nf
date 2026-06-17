/*
 * =====================================================================================
 *  modules/celloracle/perturbation/main.nf
 * -------------------------------------------------------------------------------------
 *  Process:  CO_PERTURBATION
 *  Script :  bin/co_perturbation.py
 *  Purpose:  In-silico perturbation of a single TF (default: Mafb KO) on the
 *            trained Oracle, with grid rasterisation and mass filter.
 * =====================================================================================
 */

process CO_PERTURBATION {
    tag  "perturbation_${params.perturb_tf}"
    label 'celloracle'

    input:
    path adata_pkl
    path oracle_pkl
    path links_pkl

    output:
    path "oracle_perturbed.pkl" , emit: oracle
    path "*.png"                , emit: plots, optional: true
    path "versions.yml"         , emit: versions

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/co_perturbation.py \\
        --input-adata     ${adata_pkl} \\
        --input-oracle    ${oracle_pkl} \\
        --input-links     ${links_pkl} \\
        --cluster-col     "${params.pseudotime_cluster_col}" \\
        --perturb-tf      "${params.perturb_tf}" \\
        --perturb-value   ${params.perturb_value} \\
        --n-propagation   ${params.perturb_n_propagation} \\
        --n-neighbors     ${params.perturb_n_neighbors} \\
        --sigma-corr      ${params.perturb_sigma_corr} \\
        --p-mass-n-grid   ${params.p_mass_n_grid} \\
        --p-mass-smooth   ${params.p_mass_smooth} \\
        --p-mass-min-mass ${params.p_mass_min_mass} \\
        --output-oracle   oracle_perturbed.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:    \$(python --version | sed 's/Python //')
        celloracle:\$(python -c "import celloracle; print(celloracle.__version__)")
    END_VERSIONS
    """
}
