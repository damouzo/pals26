# `modules/celloracle/perturbation/`

Single-process module `CO_PERTURBATION`. Simulates the
perturbation of a TF (default: `Mafb` knock-out) on the trained
Oracle, rasterises the result on a `n_grid × n_grid` lattice, applies
a `min_mass` filter, and pickles the perturbed Oracle for the
gradient / dev module.

See `../../../docs/steps/04_perturbation.md`.
