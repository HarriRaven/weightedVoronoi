# Weighted Euclidean tessellation (core)

Internal/core function used by
[`weighted_voronoi_domain()`](https://HarriRaven.github.io/weightedVoronoi/reference/weighted_voronoi_domain.md)
to compute a weighted Euclidean tessellation on a rasterised domain.

## Usage

``` r
weighted_voronoi(
  points_sf,
  weight_col,
  boundary = NULL,
  template_rast = NULL,
  res = NULL,
  weight_transform = function(w) w,
  method = c("argmin", "partition"),
  max_dist = NULL,
  verbose = TRUE,
  island_min_cells = 5,
  island_fill_iter = 50
)
```

## Arguments

- points_sf:

  An `sf` POINT object containing generator locations and attributes.

- weight_col:

  Character. Name of the weight column in `points_sf`.

- boundary:

  Optional `sf` polygon defining the tessellation domain. Used when
  `template_rast` is `NULL`.

- template_rast:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  template raster. Provide this instead of `boundary` + `res`.

- res:

  Numeric. Raster resolution in CRS units (e.g. metres).

- weight_transform:

  Function. Transforms weights before allocation. Must return finite,
  strictly positive values.

- method:

  Character. Allocation method; one of `"argmin"` or `"partition"`.

- max_dist:

  Optional numeric. Maximum Euclidean distance to consider (euclidean
  only).

- verbose:

  Logical. If `TRUE`, prints progress.

- island_min_cells:

  Integer. Minimum patch size used in island removal.

- island_fill_iter:

  Integer. Maximum iterations for filling reassigned cells.

## Value

A list containing polygon output (if requested), allocation raster, and
weights.
