# `modules/celloracle/pseudotime/`

Single-process module `CO_PSEUDOTIME`. Runs
`celloracle.applications.Pseudotime_calculator` on the preprocessed
AnnData, writes a pseudotime-coloured UMAP and pickles the AnnData
with the new `Pseudotime` obs column.

See `../../../docs/steps/02_pseudotime.md`.
