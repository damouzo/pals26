# pals26

> Educational Nextflow (DSL2) implementation of the
> [CellOracle](https://github.com/morris-lab/CellOracle) **scvelo-pancreas** tutorial.
> Built for the PALS 2026 course.

The pipeline takes a single-cell RNA-seq dataset (default: scvelo's `pancreas`),
infers a Gene Regulatory Network with CellOracle, performs an in-silico
perturbation (`Mafb` knock-out by default) and projects the result onto a
pseudotime-aware development field.

---

## Quick start

```bash
# Workstation (Docker)
nextflow run main.nf -profile local
```

> **Minimum hardware (laptop / `-profile local`):** 2 CPUs and 6 GB of
> free RAM. The scvelo `pancreas` dataset is small (~3.6k cells,
> ~2k HVGs) and fits comfortably inside that envelope. If you have less
> than 6 GB free, override per-run:
>
> ```bash
> nextflow run main.nf -profile local --executor.memory '4 GB'
> ```
>
> For a beefier workstation (4+ CPUs, 16+ GB) just bump the `local`
> profile block in `nextflow.config` (executor ceiling and the
> `celloracle` label).

```bash
# QMUL Apocrita HPC (SLURM + Singularity runtime)
nextflow run main.nf -profile apocrita,singularity
```

### Container on Apocrita (pre-download recommended)

Singularity's Go runtime fans out dozens of goroutines while pulling and
unpacking OCI layers. On a shared login node this can exhaust the per-user
thread limit and crash with:

```
runtime/cgo: pthread_create failed: Resource temporarily unavailable
```

The pipeline uses a local Singularity cache under `work/singularity/`
and raises `pullTimeout` to 60 min to mitigate this. If the pull still
fails, pre-download the SIF once on a compute node and point
`CELLORACLE_CONTAINER` at the local file — Nextflow will then skip the
registry entirely:

```bash
srun --cpus-per-task=4 --mem=16G --pty bash
module load nextflow
export NXF_SINGULARITY_CACHEDIR=$PWD/work/singularity
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
GOMAXPROCS=1 singularity pull \
    "$NXF_SINGULARITY_CACHEDIR/pals26-latest.sif" \
    docker://ghcr.io/damouzo/pals26:latest
export CELLORACLE_CONTAINER="$NXF_SINGULARITY_CACHEDIR/pals26-latest.sif"
nextflow run main.nf -profile apocrita,singularity
```

# Show all options
nextflow run main.nf --help
```

The pipeline is fully self-contained: the scvelo dataset and the
CellOracle base GRN are downloaded automatically on first run.

---


