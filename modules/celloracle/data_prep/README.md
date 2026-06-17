# `modules/celloracle/data_prep/`

Single-process module `CO_DATA_PREP`. Loads the scvelo dataset,
runs standard scanpy preprocessing, writes a Leiden-coloured UMAP
and pickles the AnnData for the downstream processes.

See `../../../docs/steps/01_data_prep.md` for the full explanation.
