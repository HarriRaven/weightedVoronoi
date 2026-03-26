# Uncertainty-aware weighted tessellation

Repeats weighted tessellation under stochastic perturbation of generator
weights and summarises the results as per-cell membership probabilities,
modal allocation, and entropy.

## Usage

``` r
weighted_voronoi_uncertainty(
  points_sf,
  weight_col,
  boundary_sf,
  n_sim = 100,
  weight_sd = NULL,
  distance = c("euclidean", "geodesic"),
  geodesic_engine = c("multisource", "classic"),
  res = 20,
  keep_simulations = FALSE,
  seed = NULL,
  warn_zero_entropy = TRUE,
  prepared = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- points_sf:

  An `sf` POINT object containing generator locations and attributes.

- weight_col:

  Character. Name of the weight column in `points_sf`.

- boundary_sf:

  An `sf` POLYGON/MULTIPOLYGON defining the domain.

- n_sim:

  Integer. Number of simulation runs.

- weight_sd:

  Optional numeric. Standard deviation of lognormal weight perturbation
  on the log scale. If `NULL`, no perturbation is applied and the same
  tessellation is repeated.

- distance:

  Character. One of `"euclidean"` or `"geodesic"`.

- geodesic_engine:

  Character. Geodesic engine to use when `distance = "geodesic"`.
  Defaults to `"multisource"`.

- res:

  Numeric. Raster resolution in CRS units (e.g. metres).

- keep_simulations:

  Logical. If `TRUE`, return the simulated allocation rasters as a
  stack.

- seed:

  Optional integer random seed for reproducibility.

- warn_zero_entropy:

  Logical. If `TRUE`, warn when all entropy values are zero across the
  domain.

- prepared:

  Optional prepared geodesic context created by
  [`prepare_geodesic_context()`](https://HarriRaven.github.io/weightedVoronoi/reference/prepare_geodesic_context.md)
  for repeated compatible geodesic runs.

- verbose:

  Logical. If `TRUE`, prints progress.

- ...:

  Additional arguments passed to
  [`weighted_voronoi_domain()`](https://HarriRaven.github.io/weightedVoronoi/reference/weighted_voronoi_domain.md).

## Value

A list with probability surfaces, modal allocation, entropy, and
optionally the full simulation stack.

## Details

This first implementation supports uncertainty in generator weights
only. Weights are perturbed independently across simulations using a
lognormal multiplicative model:

`w_sim = w * exp(N(0, weight_sd))`

The output includes:

- probabilities:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with one layer per generator, containing the fraction of simulations
  in which each cell was assigned to that generator.

- modal_allocation:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  giving the most probable generator for each cell.

- entropy:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  showing per-cell uncertainty. Higher values indicate less stable
  allocation across simulations.
