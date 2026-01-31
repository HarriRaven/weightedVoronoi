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
  distance = c("euclidean", "geodesic"),
  max_dist = NULL,
  island_min_cells = 5,
  island_fill_iter = 50,
  clip_to_boundary = TRUE,
  close_mask = TRUE,
  close_iters = 1,
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
```
