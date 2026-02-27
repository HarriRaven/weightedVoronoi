# Compose a resistance surface from multiple raster layers

Aligns one or more
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
layers to a common template (CRS, resolution, extent) and combines them
into a single resistance raster using a specified rule. Intended for
building `resistance_rast` inputs for geodesic tessellations.

## Usage

``` r
compose_resistance(
  ...,
  template = NULL,
  mask = NULL,
  method = c("multiply", "add", "max"),
  resample_method = c("bilinear", "near"),
  na_policy = c("propagate", "ignore")
)
```

## Arguments

- ...:

  One or more
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  layers. All layers must represent strictly positive resistance values.

- template:

  Optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  defining the target grid (CRS, resolution, extent). Defaults to the
  first layer.

- mask:

  Optional mask applied after alignment. Can be a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  or `sf` polygon.

- method:

  How to combine layers: `"multiply"` (default), `"add"`, or `"max"`.

- resample_method:

  Resampling method used when aligning layers: `"bilinear"` for
  continuous resistance, `"near"` for categorical rasters.

- na_policy:

  How to treat NA values: `"propagate"` makes any NA propagate to the
  output; `"ignore"` drops NAs (neutral element for the chosen method).

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
resistance surface on the template grid.

## Examples

``` r
if (FALSE) { # \dontrun{
R <- compose_resistance(slope_resistance, landcover_resistance, method = "multiply")
R <- add_barriers(R, rivers_sf, permeability = "semi", cost_multiplier = 20, width = 30)
out <- weighted_voronoi_domain(points_sf, "w", boundary_sf, distance = "geodesic",
                               resistance_rast = R)
} # }
```
