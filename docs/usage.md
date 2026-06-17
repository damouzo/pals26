# Usage guide

## 1. Requirements

| Tool   | Version    | Notes                                  |
|--------|------------|----------------------------------------|
| Nextflow | `>= 23.10` | `conda install -c bioconda nextflow`   |
| Docker or Singularity | latest | Apocrita uses Singularity, everything else Docker |
| A POSIX shell | bash ≥ 4 | Windows: WSL2 recommended             |

The pipeline itself downloads all data dependencies; you only need
the container runtime + Nextflow.

### Hardware

| Profile    | Min CPUs | Min free RAM | Notes                          |
|------------|----------|--------------|--------------------------------|
| `local`    | 2        | 6 GB         | Laptop / workstation           |
| `test`     | 1        | 3 GB         | Smoke run only                 |
| `apocrita` | 4        | 16 GB        | Per SLURM task                 |

If your laptop has less than 6 GB free, pass
`--executor.memory '4 GB'` (and consider downgrading
`process.withLabel: celloracle.memory` in `nextflow.config` to match).

## 2. Profiles

| Profile    | Where it runs                          | Container engine |
|------------|----------------------------------------|------------------|
| `local`    | Laptop / workstation                   | Docker           |
| `singularity` | Singularity runtime config          | Singularity      |
| `apocrita` | QMUL Apocrita SLURM cluster            | SLURM            |
| `test`     | Tiny smoke run (lower resources)       | Docker           |

```bash
nextflow run main.nf -profile local
nextflow run main.nf -profile apocrita,singularity
nextflow run main.nf -profile test
```

## 3. Common parameters

| Parameter                | Default                                | Description                                   |
|--------------------------|----------------------------------------|-----------------------------------------------|
| `--outdir`               | `./results`                            | Output directory                              |
| `--leiden_resolution`    | `0.3`                                  | Leiden clustering resolution                  |
| `--n_pca_components`     | `26`                                   | PCA components used by CellOracle             |
| `--knn_k`                | `92`                                   | k for KNN-imputation (≈0.025 × n_cells)       |
| `--grn_alpha`            | `10`                                   | Ridge α for the GRN regression                |
| `--perturb_tf`           | `Mafb`                                 | TF to perturb                                 |
| `--perturb_value`        | `0.0`                                  | Target expression of the perturbed TF         |
| `--p_mass_n_grid`        | `40`                                   | Grid resolution for the perturbation field    |
| `--p_mass_min_mass`      | `6.2`                                  | Mass cutoff for informative grid points       |
| `--celloracle_container` | `ghcr.io/<owner>/pals_celloracle:latest` | Override the published image                |
| `--help`                 | —                                      | Print the CLI help and exit                   |

Any of these can be supplied on the command line or as environment
variables (Nextflow convention: `CELLORACLE_CONTAINER` etc.).

## 4. Skipping / re-running steps

Every step can be disabled with `--run_<step> false`. When a step is
disabled, the pipeline expects the corresponding output file in
`--outdir/<step_dir>/`. This lets you iterate on, e.g., the
perturbation step without re-running the slow GRN inference:

```bash
# Re-run only the perturbation + gradient steps
nextflow run main.nf -profile local \
    --run_preprocess  false \
    --run_pseudotime  false \
    --run_grn         false
```

## 5. Apocrita (SLURM) tips

* The repo does not hardcode a Slurm account. If your site needs one,
  add it in a local config override or pass site-specific `clusterOptions`.
* Singularity images are cached under `work/singularity/` by default.
  If you want a different path, set `NXF_SINGULARITY_CACHEDIR` before
  launching.
* For very large runs, increase `process.cpus` / `process.memory` in
  `conf/base.config` or override per-process via `withName:` in
  `conf/modules.config`.

## 6. Development tips

* Each Python script is fully isolated under `bin/` and uses only
  `argparse`. You can run them directly for debugging:

  ```bash
  python bin/co_data_prep.py --help
  ```
* To test a single step end-to-end, comment out the `if (...)` blocks
  in `subworkflows/local/celloracle.nf` and re-run with `-resume`.
* `versions.yml` files are emitted by every process so you can audit
  exactly which package versions were used.

## 7. Help & support

* Run `nextflow run main.nf --help` for a quick reminder.
* Open an issue on GitHub: <https://github.com/BCI-KRP/pals_celloracle/issues>
* Nextflow documentation: <https://www.nextflow.io/docs/latest/index.html>
