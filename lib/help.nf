/*
 * =====================================================================================
 *  lib/help.nf
 * -------------------------------------------------------------------------------------
 *  Shared CLI helpers used by main.nf. Kept out of the entry-point so
 *  main.nf can read like a one-liner: include + workflow { CELLORACLE() }.
 * =====================================================================================
 */

def helpMessage() {
    log.info"""
    =============================================================
     pals_celloracle v${workflow.manifest.version} — CellOracle tutorial, Nextflow-style
    =============================================================

     Required arguments:
       (none) — the scvelo pancreas dataset and CellOracle base GRN
                are downloaded automatically inside the pipeline.

     Pipeline switches:
       --run_preprocess        Run the preprocessing + Leiden + UMAP step     [default: true]
       --run_pseudotime        Run the diffusion-pseudotime step             [default: true]
       --run_grn               Run the CellOracle GRN-inference step         [default: true]
       --run_perturbation      Run the in-silico perturbation step            [default: true]
       --run_gradient          Run the gradient + development-module step     [default: true]

     Algorithm tuning (defaults match the upstream tutorial; set to "auto" to derive from data):
       --leiden_resolution     Leiden resolution                             [default: 0.3]
       --n_pca_components      PCA components ("auto" => knee point)         [default: auto]
       --knn_k                 KNN k ("auto" => 2.5% of cells)               [default: auto]
       --grn_alpha             Ridge alpha for GRN regression                [default: 10]
       --perturb_tf            TF to knock-out                               [default: Mafb]
       --perturb_value         Target expression of the perturbed TF         [default: 0.0]
       --p_mass_n_grid         Grid resolution for the perturbation field    [default: 40]
       --p_mass_min_mass       Mass cutoff to keep informative grid points   [default: 6.2]

     I/O & runtime:
       --outdir                Output directory                              [default: ./results]
       --publish_dir_mode      'copy' | 'symlink' | 'link'                   [default: copy]
       --celloracle_container  Override the container image                  [default: ghcr.io/damouzo/pals_celloracle:latest]

     Profiles:
       local        Run on a single workstation (Docker).
      apocrita     Run on QMUL's Apocrita cluster (SLURM only).
      singularity  Singularity runtime settings and cache handling.
       test         Tiny smoke run with reduced resources.

     Help & documentation:
       --help                 Show this help message and exit
       docs/usage.md          Full user-guide and per-process explanations
       docs/steps/            Per-process step-by-step notes for students
       https://www.nextflow.io/docs/latest/index.html  — Nextflow open-source docs
    """.stripIndent()
}

def logRunSummary() {
    log.info """
    =============================================================
     pals_celloracle v${workflow.manifest.version}
    =============================================================
     Output directory   : ${params.outdir}
     Container          : ${params.celloracle_container}
     Steps enabled      : preprocess=${params.run_preprocess}, pseudotime=${params.run_pseudotime},
                          grn=${params.run_grn}, perturbation=${params.run_perturbation},
                          gradient=${params.run_gradient}
     GRN                : n_pca=${params.n_pca_components}, k=${params.knn_k}, alpha=${params.grn_alpha}
     Perturbation       : TF=${params.perturb_tf} -> ${params.perturb_value},
                          n_propagation=${params.perturb_n_propagation},
                          n_grid=${params.p_mass_n_grid}, min_mass=${params.p_mass_min_mass}
     Profile            : ${workflow.profile}
    =============================================================
    """
}

def containerPullHint() {
    def resolved = (params.celloracle_container ?: '').toString()
    if (resolved.startsWith('docker://') || resolved.startsWith('ghcr.io/') || resolved.startsWith('docker.io/')) {
        log.info """
        =============================================================
         NOTE — Container will be pulled from the registry
        =============================================================
         Image : ${resolved}
         To use a local image instead:
           nextflow run main.nf --celloracle_container /path/to/image.sif
         or via the generic Nextflow flag:
           nextflow run main.nf -with-container /path/to/image.sif
        =============================================================
        """.stripIndent()
    }
}
