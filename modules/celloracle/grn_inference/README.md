# `modules/celloracle/grn_inference/`

Single-process module `CO_GRN_INFERENCE`. Builds a `co.Oracle`,
imports the mouse scATAC base GRN, performs KNN-imputation and fits
a per-cluster ridge regression. The `Oracle` + `Links` objects are
pickled for the perturbation step.

This module forwards `"auto"` for `--n-pca` and `--knn-k` so that
`bin/co_grn_inference.py` derives them from the data
(knee point / 2.5% of cells).

See `../../../docs/steps/03_grn_inference.md`.
