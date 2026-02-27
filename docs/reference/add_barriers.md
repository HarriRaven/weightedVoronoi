# Add barriers to a resistance surface

Modifies a resistance raster by applying semi-permeable or impermeable
barriers provided as a raster mask (values \> 0 treated as barrier) or
as vector features (sf / SpatVector) rasterised onto the resistance
grid.

## Usage

``` r
add_barriers(
  resistance,
  barriers,
  permeability = c("semi", "impermeable", "permeable"),
  cost_multiplier = 10,
  width = 0
)
```

## Arguments

- resistance:

  terra::SpatRaster of strictly positive movement resistance.

- barriers:

  A terra::SpatRaster mask (values \> 0 treated as barrier), or an
  sf/SpatVector LINESTRING/POLYGON object.

- permeability:

  One of "semi", "impermeable", "permeable".

- cost_multiplier:

  Numeric \> 0. Multiplier applied where barrier present (for "semi" and
  "permeable").

- width:

  Buffer distance (CRS units) applied to vector barriers before
  rasterising.

## Value

A terra::SpatRaster resistance surface with barrier effects applied.
