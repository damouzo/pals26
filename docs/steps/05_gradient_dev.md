# Step 5 — `CO_GRADIENT_DEV` — Gradient & development module

**Process file:** `modules/celloracle_steps.nf`
**Script:** `bin/co_gradient_dev.py`

## What it does

The final step compares the **simulated** perturbation flow with a
**reference** differentiation flow (built from the diffusion
pseudotime) to highlight the regions where the perturbation pushes
cells *away* from their normal developmental trajectory.

1. Loads the perturbed Oracle and injects the `Pseudotime` column
   saved in step 2.
2. Builds a `Gradient_calculator` on the **same grid** as the
   perturbation (`n_grid=40`, `n_neighbors=200`).
3. `gradient.calculate_p_mass` + `gradient.calculate_mass_filter`
   use the same `min_mass=6.2` cutoff, so both fields cover the
   same set of grid points.
4. `gradient.transfer_data_into_grid(method='polynomial', n_poly=3)`
   fits a 3rd-degree polynomial of the pseudotime on the grid.
5. `gradient.calculate_gradient()` yields the reference
   differentiation vectors.
6. `Oracle_development_module` loads both the gradient and the
   perturbation, then computes the **inner product** of the two
   vector fields:
   - **positive** inner product → perturbation *follows* normal
     differentiation
   - **negative** inner product → perturbation *diverts* cells from
     the normal trajectory
7. `dev.calculate_digitized_ip(n_bins=10)` bins the inner-product
   values into 10 quantile bins per cell.
8. Renders the headline plot (`pertubation_score.png`) that overlays
   the inner-product colour field with the perturbation flow arrows.
9. Writes the per-cell digitised scores to `dev_scores.tsv`.

## Inputs

| Channel | File | Description |
|---------|------|-------------|
| `input-adata`  | `adata_with_pseudotime.pkl` | Supplies the `Pseudotime` column |
| `input-oracle` | `oracle_perturbed.pkl`     | Oracle after the simulation |

## Outputs

| File | Description |
|------|-------------|
| `pertubation_score.png` | Inner-product field + perturbation flow |
| `dev_scores.tsv`        | Per-cell digitised inner-product scores  |
| `versions.yml`          | Tool versions                             |

## Knobs

| Flag               | Default | Effect                              |
|--------------------|---------|-------------------------------------|
| `--pseudotime-key` | `Pseudotime` | `obs` column carrying pseudotime |
| `--n-grid`         | 40      | Must match `--p_mass_n_grid` from step 4 |
| `--n-poly`         | 3       | Polynomial degree for grid transfer |
| `--n-bins`         | 10      | Number of digitisation bins         |

## Reading the final plot

* **Bright yellow** cells: large *positive* inner product — the
  perturbation moves cells in the same direction as normal
  development.
* **Dark purple** cells: large *negative* inner product — the
  perturbation pushes cells in the opposite direction of normal
  development.
* **Arrows**: the 2-D flow field produced by the in-silico
  perturbation; their size is scaled by `--scale=30`.

For the Mafb knock-out in the pancreas dataset, the brightest
"diverted" cells are expected to sit in the late-β / late-α branch,
matching the published CellOracle result.
