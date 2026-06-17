/*
 * =====================================================================================
 *  modules/celloracle/gradient_dev/main.nf
 * -------------------------------------------------------------------------------------
 *  Process:  CO_GRADIENT_DEV
 *  Script :  bin/co_gradient_dev.py
 *  Purpose:  Build the reference differentiation gradient and the
 *            Oracle_development_module, render the headline inner-product plot.
 * =====================================================================================
 */

process CO_GRADIENT_DEV {
    tag  "gradient_dev"
    label 'celloracle'

    input:
    path adata_pkl
    path oracle_pkl

    output:
    path "pertubation_score.png" , emit: plot
    path "dev_scores.tsv"        , emit: scores
    path "versions.yml"          , emit: versions

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/co_gradient_dev.py \\
        --input-adata    ${adata_pkl} \\
        --input-oracle   ${oracle_pkl} \\
        --pseudotime-key "Pseudotime" \\
        --n-grid         ${params.p_mass_n_grid} \\
        --n-neighbors    ${params.perturb_n_neighbors} \\
        --smooth         ${params.p_mass_smooth} \\
        --min-mass       ${params.p_mass_min_mass} \\
        --n-poly         ${params.gradient_n_poly} \\
        --n-bins         ${params.dev_n_bins} \\
        --output-plot    pertubation_score.png \\
        --output-scores  dev_scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:     \$(python --version | sed 's/Python //')
        celloracle: \$(python -c "import celloracle; print(celloracle.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
