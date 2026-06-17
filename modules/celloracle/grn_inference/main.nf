/*
 * =====================================================================================
 *  modules/celloracle/grn_inference/main.nf
 * -------------------------------------------------------------------------------------
 *  Process:  CO_GRN_INFERENCE
 *  Script :  bin/co_grn_inference.py
 *  Purpose:  Train a CellOracle Oracle + Links object on the cluster-level GRN.
 *            Supports "auto" for --n-pca and --knn-k so the script derives them
 *            from the data shape (knee point of explained variance; 2.5% of cells).
 * =====================================================================================
 */

process CO_GRN_INFERENCE {
    tag  "grn"
    label 'celloracle'

    input:
    path adata_pkl

    output:
    path "oracle_object.pkl" , emit: oracle
    path "links_object.pkl"  , emit: links
    path "grn_score.png"     , emit: plot
    path "grn_*.png"         , emit: extra_plots, optional: true
    path "versions.yml"      , emit: versions

    script:
    // Forward an "auto" sentinel to Python when the user did not pin a value.
    def pca_dim = (params.n_pca_components ?: "auto").toString()
    def knn_k   = (params.knn_k           ?: "auto").toString()

    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/co_grn_inference.py \\
        --input-adata   ${adata_pkl} \\
        --cluster-col   "${params.pseudotime_cluster_col}" \\
        --n-pca         "${pca_dim}" \\
        --knn-k         "${knn_k}" \\
        --knn-n-jobs    ${task.cpus} \\
        --grn-alpha     ${params.grn_alpha} \\
        --filter-p      ${params.grn_filter_p} \\
        --filter-max    ${params.grn_filter_max_links} \\
        --output-oracle oracle_object.pkl \\
        --output-links  links_object.pkl \\
        --output-plot   grn_score.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:    \$(python --version | sed 's/Python //')
        celloracle:\$(python -c "import celloracle; print(celloracle.__version__)")
        scanpy:    \$(python -c "import scanpy;    print(scanpy.__version__)")
    END_VERSIONS
    """
}
