# =====================================================================================
#  pals_celloracle — reproducible environment
# -------------------------------------------------------------------------------------
#  Python 3.8 is pinned because CellOracle (latest published release) is built and
#  tested against 3.8; later Python versions break the bundled C extensions.
#  All heavy install steps live under a single RUN to minimise layers and image size.
# =====================================================================================
FROM condaforge/miniforge3:latest AS build

LABEL org.opencontainers.image.title="pals_celloracle" \
      org.opencontainers.image.description="Reproducible Nextflow environment for the CellOracle scvelo-pancreas tutorial" \
      org.opencontainers.image.source="https://github.com/BCI-KRP/pals_celloracle" \
      org.opencontainers.image.licenses="MIT"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# -------------------------------------------------------------------------------------
#  System libraries required by leidenalg, louvain, scvelo, h5py, matplotlib
# -------------------------------------------------------------------------------------
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        gfortran \
        git \
        wget \
        ca-certificates \
        libhdf5-dev \
        libxml2-dev \
        libxslt-dev \
        libpng-dev \
        libjpeg-dev \
        libgomp1 \
 && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------------------------
#  Conda environment with the full scientific Python + scRNA-seq stack
#  Pre-installing heavy libraries via conda prevents compiling from source via pip.
# -------------------------------------------------------------------------------------
RUN mamba create -n celloracle -y -c conda-forge -c bioconda \
        python=3.8 \
        numpy=1.23 \
        pandas=1.5 \
        scipy=1.10 \
        scikit-learn=1.1 \
        matplotlib=3.5 \
        seaborn=0.12 \
        h5py=3.7 \
        networkx=2.8 \
        joblib=1.2 \
        tqdm=4.64 \
        numba=0.56 \
        umap-learn=0.5 \
        statsmodels=0.13 \
        leidenalg=0.9 \
        igraph=0.10 \
        louvain=0.7 \
        scanpy=1.9.3 \
        scvelo=0.2.5 \
        anndata=0.9.2 \
 && mamba clean -afy

# Make the conda env the default Python
ENV PATH=/opt/conda/envs/celloracle/bin:$PATH

# -------------------------------------------------------------------------------------
#  pip-only packages — CellOracle via PyPI to avoid broken GitHub Git Checkouts
# -------------------------------------------------------------------------------------
RUN pip install --no-cache-dir -U setuptools wheel pip \
 && pip install --no-cache-dir "cython<3" \
 && pip install --no-cache-dir \
        scrublet==0.2.3 \
        harmony-pytorch==0.1.6 \
        fa2==0.3.5 \
        celloracle==0.18.0

# CellOracle downloads a base GRN the first time it is requested; cache it at build
# time so the pipeline does not need network access on first run.
RUN python -c "import celloracle as co; co.data.load_mouse_scATAC_atlas_base_GRN()"

# -------------------------------------------------------------------------------------
#  Final stage: slim runtime image
FROM condaforge/miniforge3:latest

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# Copy the prebuilt env (faster than re-installing in the runtime image)
COPY --from=build /opt/conda /opt/conda

# Minimal system libraries for the runtime image
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        libgomp1 \
        libpng16-16 \
        libjpeg-turbo-progs \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/envs/celloracle/bin:$PATH

# Smoke test: import the heavy libraries to fail fast at build time
RUN python -c "import scanpy, scvelo, celloracle, leidenalg; print('celloracle', celloracle.__version__); print('scanpy', scanpy.__version__); print('scvelo', scvelo.__version__)"

WORKDIR /work
CMD ["/bin/bash"]