# Step 4 — `CO_PERTURBATION` — In-silico perturbation

**Process file:** `modules/celloracle_steps.nf`
**Script:** `bin/co_perturbation.py`

## What it does

Takes the trained Oracle + Links and simulates the effect of
perturbing a single transcription factor (default: knock-out of
**Mafb**):

1. `Oracle.get_cluster_specific_TFdict_from_Links(links_object=links)`
   — restricts the perturbation model to TFs that actually have an
   outgoing edge in each cluster (avoids testing TFs that the GRN
   says do nothing).
2. `Oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)`
   — fits a new regression model with the same `α` as the inference
   step, this time restricted to the cluster-specific TF set.
3. `Oracle.simulate_shift(perturb_condition={'Mafb': 0.0},
   n_propagation=3)` — propagates the perturbation through three
   GRN iterations, yielding a delta-expression vector per cell.
4. `Oracle.estimate_transition_prob(n_neighbors=200, knn_random=True,
   sampled_fraction=1)` — builds a Markov transition matrix on the
   embedding to translate the expression shift into a 2-D flow.
5. `Oracle.calculate_embedding_shift(sigma_corr=0.05)` — converts
   expression deltas to UMAP-space vectors.
6. `Oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)`
   — rasterises the simulated cells onto a 40×40 grid.
7. `Oracle.calculate_mass_filter(min_mass=6.2, plot=False)` — masks
   out grid points that contain fewer than 6.2 cells worth of
   information (avoids noisy, isolated bins).
8. Pickles the perturbed Oracle.

## Inputs

| Channel | File | Description |
|---------|------|-------------|
| `input-adata`  | `adata_with_pseudotime.pkl` | Carries the Leiden labels (used by `fit_GRN_for_simulation`) |
| `input-oracle` | `oracle_object.pkl`        | Trained Oracle from step 3 |
| `input-links`  | `links_object.pkl`         | Per-cluster TF list (step 3) |

## Outputs

| File | Description |
|------|-------------|
| `oracle_perturbed.pkl` | Oracle with simulated shift + flow field |
| `versions.yml`         | Tool versions                              |

## Knobs

| Flag                  | Default | Effect                                          |
|-----------------------|---------|-------------------------------------------------|
| `--perturb-tf`        | `Mafb`  | TF to perturb                                   |
| `--perturb-value`     | `0.0`   | New expression of the TF (0.0 ≈ full knock-out) |
| `--n-propagation`     | 3       | GRN propagation rounds                          |
| `--n-neighbors`       | 200     | Neighbours for the KNN flow                     |
| `--sigma-corr`        | 0.05    | Correlation-σ for embedding shift               |
| `--p-mass-n-grid`     | 40      | Grid resolution                                 |
| `--p-mass-smooth`     | 0.8     | Smoothing factor for `calculate_p_mass`         |
| `--p-mass-min-mass`   | 6.2     | Mass cutoff (filters out sparse bins)           |

## Sanity check

The perturbed Oracle's `delta_embedding` attribute holds the
per-cell UMAP displacement vectors. The next process (`CO_GRADIENT_DEV`)
loads this directly from the pickle to compute inner products with
the reference differentiation flow.
