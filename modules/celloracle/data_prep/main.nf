/*
 * =====================================================================================
 *  modules/celloracle/data_prep/main.nf
 * -------------------------------------------------------------------------------------
 *  Process:  CO_DATA_PREP
 *  Script :  bin/co_data_prep.py
 *  Purpose:  Load the scvelo pancreas dataset (auto-downloaded by scv.datasets),
 *            normalise / log-transform / select HVGs, run PCA + Leiden + UMAP,
 *            pickle the AnnData for downstream processes.
 * =====================================================================================
 */

process CO_DATA_PREP {
    tag  "${dataset}"
    label 'celloracle'

    input:
    val dataset

    output:
    path "adata_preprocessed.pkl" , emit: adata
    path "qc_umap.png"             , emit: umap
    path "celltype_umap.png"       , emit: celltype_umap, optional: true
    path "versions.yml"            , emit: versions

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    python ${projectDir}/bin/co_data_prep.py \\
        --dataset           "${dataset}" \\
        --leiden-resolution ${params.leiden_resolution} \\
        --filter-min-cells  ${params.filter_min_cells} \\
        --output-adata      adata_preprocessed.pkl \\
        --output-umap       qc_umap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python:    \$(python --version | sed 's/Python //')
        scanpy:    \$(python -c "import scanpy;    print(scanpy.__version__)")
        scvelo:    \$(python -c "import scvelo;    print(scvelo.__version__)")
        anndata:   \$(python -c "import anndata;   print(anndata.__version__)")
        leidenalg: \$(python -c "import leidenalg; print(leidenalg.version)")
    END_VERSIONS
    """
}
