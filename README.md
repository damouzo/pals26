# pals26

Nextflow of the CellOracle tutorial for the PALS 2026 course. It takes a scRNA-seq dataset, infers a Gene Regulatory Network with CellOracle, runs an in-silico knock-out (Mafb by default) and projects the result onto a pseudotime-aware development field.

## Requirements

- [Git](https://git-scm.com/install/)
- Java 17 or newer
- [Nextflow](https://docs.seqera.io/nextflow/install) (>= 23.10.0)
- A container runtime: [Docker](https://docs.docker.com/engine/install/) for your laptop, or Singularity/Apptainer if you are running on an HPC

## Install

### Nextflow

Nextflow needs Java 17. If you don't have it yet, the easiest way is via [SDKMAN!](https://sdkman.io/):

```bash
curl -s https://get.sdkman.io | bash   # then open a new terminal
sdk install java 17.0.10-tem
```

Then install Nextflow itself:

```bash
curl -s https://get.nextflow.io | bash
mkdir -p $HOME/.local/bin
mv nextflow $HOME/.local/bin/
export PATH="$PATH:$HOME/.local/bin"   # add this to your ~/.bashrc or ~/.zshrc
```

Check that everything is wired up:

```bash
java -version
nextflow -version
```

See the [official install guide](https://docs.seqera.io/nextflow/install) for other options (Conda, standalone, Windows/WSL, etc.).

## Run

### On your laptop (Docker)

```bash
nextflow run main.nf -profile local
```

The `local` profile in `nextflow.config` defaults to 2 CPUs and 6 GB of RAM, which is enough for the scvelo `pancreas` dataset. If your machine is beefier, bump the `executor` and `celloracle` label blocks in `nextflow.config`. If you are tight on RAM, lower it per-run:

```bash
nextflow run main.nf -profile local --executor.memory '4 GB'
```

### On a SLURM cluster (Singularity/Apptainer)

```bash
nextflow run main.nf -profile cluster,singularity
```

The `cluster` profile uses SLURM with 4 CPUs, 32 GB RAM and a 4 h walltime per task. Tweak `nextflow.config` to match your scheduler and quotas.

## Data

The scvelo `pancreas` dataset and the CellOracle base GRN are downloaded automatically on the first run, so you don't need to fetch anything by hand.

## Help

```bash
nextflow run main.nf --help
```

More detail in [`docs/usage.md`](docs/usage.md) and the per-step notes under [`docs/steps/`](docs/steps/).
