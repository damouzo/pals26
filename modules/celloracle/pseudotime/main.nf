/*
 * =====================================================================================
 *  modules/celloracle/pseudotime/main.nf
 * -------------------------------------------------------------------------------------
 *  Process:  CO_PSEUDOTIME
 *  Script :  bin/co_pseudotime.py
 *  Purpose:  Diffusion-pseudotime via CellOracle's Pseudotime_calculator.
 * =====================================================================================
 */

process CO_PSEUDOTIME {
    tag  "pseudotime"
    label 'celloracle'

    input:
    path adata_pkl

    output:
    path "adata_with_pseudotime.pkl" , emit: adata
    path "pt.png"                    , emit: plot
    path "versions.yml"              , emit: versions

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/co_pseudotime.py \\
        --input-adata   ${adata_pkl} \\
        --cluster-col   "${params.pseudotime_cluster_col}" \\
        --lineage-name  "${params.pseudotime_lineage}" \\
        --output-adata  adata_with_pseudotime.pkl \\
        --output-plot   pt.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:    \$(python --version | sed 's/Python //')
        scanpy:    \$(python -c "import scanpy;    print(scanpy.__version__)")
        celloracle:\$(python -c "import celloracle; print(celloracle.__version__)")
    END_VERSIONS
    """
}
