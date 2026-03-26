# Weighted tessellation in a constrained polygon domain

Creates a complete, connected tessellation of a polygonal domain using
either weighted Euclidean distance or weighted geodesic
(domain-constrained) distance. Weights are supplied as an attribute of
generator points and can be transformed by a user-defined function prior
to allocation.

## Usage

``` r
weighted_voronoi_domain(
  points_sf,
  weight_col,
  boundary_sf,
  res = 20,
  weight_transform = function(w) w,
  weight_model = c("multiplicative", "power", "additive"),
  weight_power = 1,
  distance = c("euclidean", "geodesic"),
  max_dist = NULL,
  island_min_cells = 5,
  island_fill_iter = 50,
  clip_to_boundary = TRUE,
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
  geodesic_engine = c("classic", "multisource"),
  prepared = NULL,
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

- distance:

  Character. One of `"euclidean"` or `"geodesic"`.

- max_dist:

  Optional numeric. Maximum Euclidean distance to consider (euclidean
  only).

- island_min_cells:

  Integer. Minimum patch size used in island removal.

- island_fill_iter:

  Integer. Maximum iterations for filling reassigned cells.

- clip_to_boundary:

  Logical. If `TRUE`, polygon output is intersected with the input
  boundary for exact edge matching (euclidean only).

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

- anisotropy:

  Character. Directional cost model for geodesic distance.

  "none"

  :   Standard isotropic geodesic distance (default).

  "terrain"

  :   Direction-dependent movement based on terrain slope (DEM
      required).

- uphill_factor:

  Numeric \> 0. Multiplier controlling additional cost of uphill
  movement when `anisotropy = "terrain"`. Values \> 1 penalise uphill
  movement more strongly.

- downhill_factor:

  Numeric \> 0. Multiplier controlling ease of downhill movement when
  `anisotropy = "terrain"`. Values \> 1 make downhill travel easier.

- geodesic_engine:

  Character. Geodesic allocation engine to use when
  `distance = "geodesic"`.

  "classic"

  :   Per-generator accumulated-cost allocation. Supports all current
      geodesic modes and weight models.

  "multisource"

  :   Single-pass multisource allocation. Currently supported only for
      `weight_model = "additive"` and `anisotropy = "none"`.

- prepared:

  Optional prepared geodesic context created by
  [`prepare_geodesic_context()`](https://HarriRaven.github.io/weightedVoronoi/reference/prepare_geodesic_context.md)
  for repeated compatible geodesic runs.

- verbose:

  Logical. If `TRUE`, prints progress.

## Value

A list with elements including:

- polygons:

  An `sf` object with one polygon per generator.

- allocation:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  assigning each cell to a generator.

- summary:

  A generator-level summary table.

- diagnostics:

  A list of diagnostic metrics and settings.

For geodesic allocation, `geodesic_engine = "classic"` computes one
accumulated-cost surface per generator and assigns each raster cell to
the minimum effective cost. This is the reference implementation and
supports all current geodesic modes.

`geodesic_engine = "multisource"` provides a scalable alternative for
additive-weight isotropic geodesic tessellations. It uses a single
multisource shortest-path propagation and is currently available only
when `weight_model = "additive"` and `anisotropy = "none"`.

## Details

When `distance = "geodesic"`, distances are computed as shortest paths
constrained to the spatial domain. If `dem_rast` is supplied and
`use_tobler = TRUE`, movement cost between adjacent raster cells is
modified using Tobler's hiking function, such that steeper slopes
increase effective distance. This allows elevation or resistance
surfaces to influence spatial allocation while preserving a complete
tessellation. When `distance = "geodesic"` and `anisotropy = "terrain"`,
movement costs are computed using a direction-dependent extension of a
Tobler-like hiking function.

Movement between raster cells becomes asymmetric: uphill and downhill
transitions have different costs. This results in anisotropic
(direction-dependent) geodesic tessellations.

Currently, anisotropic terrain mode:

- requires a `dem_rast` input

- does not combine with a user-supplied `resistance_rast`

- uses 8-directional neighbourhood transitions

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
crs_use <- 32636
boundary_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(rbind(
    c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
  )))),
  crs = crs_use
)
points_sf <- st_sf(
  population = c(50, 200, 1000),
  geometry = st_sfc(st_point(c(200,200)), st_point(c(800,250)), st_point(c(500,500))),
  crs = crs_use
)
out <- weighted_voronoi_domain(points_sf, "population", boundary_sf,
  res = 20, weight_transform = log10, distance = "euclidean", verbose = FALSE
)
} # }
if (FALSE) { # \dontrun{
library(sf)
library(terra)

crs_use <- "EPSG:3857"

boundary_sf <- st_sf(
  id = 1,
  geometry = st_sfc(
    st_polygon(list(rbind(
      c(0, 0), c(1000, 0), c(1000, 1000),
      c(0, 1000), c(0, 0)
    ))),
    crs = crs_use
  )
)

points_sf <- st_sf(
  population = c(1, 1),
  geometry = st_sfc(
    st_point(c(200, 500)),
    st_point(c(800, 500)),
    crs = crs_use
  )
)

dem_rast <- rast(
  ext = ext(0, 1000, 0, 1000),
  resolution = 50,
  crs = crs_use
)

xy <- crds(dem_rast, df = TRUE)
values(dem_rast) <- xy$x * 20

out <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  distance = "geodesic",
  dem_rast = dem_rast,
  anisotropy = "terrain",
  uphill_factor = 3,
  downhill_factor = 1.2
)
} # }
```
