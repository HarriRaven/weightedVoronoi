# Temporal weighted tessellation

Runs weighted tessellation across a sequence of time-specific point
datasets and returns a stack of allocation rasters, with optional
polygons and summaries per time step.

## Usage

``` r
weighted_voronoi_time(
  points_list,
  weight_col,
  boundary_sf,
  time_index = NULL,
  distance = c("euclidean", "geodesic"),
  geodesic_engine = c("multisource", "classic"),
  res = 20,
  resistance_list = NULL,
  dem_list = NULL,
  keep_polygons = FALSE,
  keep_summaries = TRUE,
  prepared = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- points_list:

  A non-empty list of `sf` POINT objects, one per time step.

- weight_col:

  Character. Name of the weight column present in each element of
  `points_list`.

- boundary_sf:

  An `sf` POLYGON/MULTIPOLYGON defining the domain.

- time_index:

  Optional character vector of time labels. Defaults to
  `names(points_list)` if present, otherwise `"t1"`, `"t2"`, etc.

- distance:

  Character. One of `"euclidean"` or `"geodesic"`.

- geodesic_engine:

  Character. Geodesic engine to use when `distance = "geodesic"`.

- res:

  Numeric. Raster resolution in CRS units (e.g. metres).

- resistance_list:

  Optional list of resistance rasters, either length 1 (reused for all
  times) or the same length as `points_list`.

- dem_list:

  Optional list of DEM rasters, either length 1 (reused for all times)
  or the same length as `points_list`.

- keep_polygons:

  Logical. If `TRUE`, return polygons for each time step.

- keep_summaries:

  Logical. If `TRUE`, return summaries for each time step.

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

A list containing:

- allocations:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with one allocation layer per time step.

- time_index:

  Character vector of time labels.

- change_map_first_last:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  indicating whether allocation changed between the first and last time
  step (`1` = changed, `0` = unchanged).

- persistence:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  indicating whether each cell retained the same allocation across all
  time steps (`1` = persistent, `0` = changed at least once).

- polygons:

  Optional list of `sf` polygon outputs by time.

- summaries:

  Optional list of summary tables by time.

## Details

This first implementation assumes a static boundary and runs each time
step independently. Time-varying weights, point locations, and
resistance/DEM surfaces are supported by supplying separate inputs for
each time step.
