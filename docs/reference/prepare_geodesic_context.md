# Prepare a geodesic context for repeated runs

Precomputes the domain mask, aligned resistance handling, transition
object, and multisource graph representation (when applicable) for
repeated geodesic tessellation workflows.

## Usage

``` r
prepare_geodesic_context(
  boundary_sf,
  res = 20,
  close_mask = TRUE,
  close_iters = 1,
  resistance_rast = NULL,
  dem_rast = NULL,
  use_tobler = TRUE,
  tobler_v0_kmh = 6,
  tobler_a = 3.5,
  tobler_b = 0.05,
  min_speed_kmh = 0.25,
  anisotropy = c("none", "terrain"),
  uphill_factor = 1,
  downhill_factor = 1,
  geodesic_engine = c("classic", "multisource")
)
```

## Arguments

- boundary_sf:

  An `sf` POLYGON/MULTIPOLYGON defining the domain.

- res:

  Numeric. Raster resolution in CRS units (e.g. metres).

- close_mask:

  Logical. If `TRUE`, applies a morphological closing to the raster
  mask.

- close_iters:

  Integer. Number of closing iterations.

- resistance_rast:

  Optional SpatRaster giving movement resistance (\>0).

- dem_rast:

  Optional SpatRaster providing elevation or resistance surface.

- use_tobler:

  Logical; if `TRUE`, apply Tobler's hiking function to convert slope
  into isotropic movement cost.

- tobler_v0_kmh:

  Base walking speed on flat terrain (km/h).

- tobler_a:

  Tobler exponential slope coefficient.

- tobler_b:

  Tobler slope multiplier.

- min_speed_kmh:

  Minimum allowed speed to avoid infinite costs.

- anisotropy:

  Character. One of `"none"` or `"terrain"`.

- uphill_factor:

  Numeric \> 0. Additional uphill movement penalty when
  `anisotropy = "terrain"`.

- downhill_factor:

  Numeric \> 0. Relative ease of downhill movement when
  `anisotropy = "terrain"`.

- geodesic_engine:

  Character. One of `"classic"` or `"multisource"`.

## Value

A prepared geodesic context object for repeated geodesic allocation.
