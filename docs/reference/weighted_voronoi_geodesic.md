# Weighted geodesic tessellation (core)

Computes a weighted tessellation using domain-constrained (geodesic)
distances. Distances are calculated as shortest-path distances through a
rasterised domain mask.

## Usage

``` r
weighted_voronoi_geodesic(
  points_sf,
  weight_col,
  boundary_sf,
  res = 20,
  weight_transform = function(w) w,
  weight_model = c("multiplicative", "power", "additive"),
  weight_power = 1,
  close_mask = TRUE,
  close_iters = 1,
  resistance_rast = NULL,
  dem_rast = NULL,
  use_tobler = TRUE,
  tobler_v0_kmh = 6,
  tobler_a = 3.5,
  tobler_b = 0.05,
  min_speed_kmh = 0.25,
  island_min_cells = 5,
  island_fill_iter = 50,
  verbose = TRUE
)
```

## Arguments

- points_sf:

  An `sf` POINT object containing generator locations and attributes.

- weight_col:

  Character. Name of the weight column in `points_sf`.

- boundary_sf:

  An `sf` POLYGON/MULTIPOLYGON defining the domain.

- res:

  Numeric. Raster resolution in CRS units (e.g. metres).

- weight_transform:

  Function. Transforms weights before allocation. Must return finite,
  strictly positive values.

- weight_model:

  Character. One of "multiplicative", "power", or "additive". Controls
  how distances and weights combine into effective cost.

- weight_power:

  Numeric \> 0. Only used when weight_model = "power". Controls the
  distance exponent.

- close_mask:

  Logical. If `TRUE`, applies a morphological closing to the raster mask
  (geodesic only).

- close_iters:

  Integer. Number of closing iterations (geodesic only).

- resistance_rast:

  Optional SpatRaster giving movement resistance (\>0). Overrides
  dem_rast/Tobler when provided.

- dem_rast:

  Optional SpatRaster providing elevation or resistance surface. Must
  align with the tessellation domain and resolution.

- use_tobler:

  Logical; if TRUE, apply Tobler's hiking function to convert slope into
  isotropic movement cost.

- tobler_v0_kmh:

  Base walking speed on flat terrain (km/h).

- tobler_a:

  Tobler exponential slope coefficient (default -3.5).

- tobler_b:

  Tobler slope multiplier (default 0.05).

- min_speed_kmh:

  Minimum allowed speed to avoid infinite costs.

- island_min_cells:

  Integer. Minimum patch size used in island removal.

- island_fill_iter:

  Integer. Maximum iterations for filling reassigned cells.

- verbose:

  Logical. If `TRUE`, prints progress.

## Value

A list containing polygon output, allocation raster, and weights.
