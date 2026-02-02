# Weighted spatial tessellations in constrained domains

``` r
library(weightedVoronoi)
```

## Introduction

Spatial allocation around point locations is a recurring requirement in
ecology and social–ecological research, including the delineation of
settlement influence zones, species territories, and management areas.
Voronoi tessellations provide a simple and transparent framework for
dividing space among generators, but standard implementations assume
uniform influence and unconstrained Euclidean distance.

The weightedVoronoi package extends this framework by allowing (i)
generator-specific weights and (ii) tessellations constrained to
arbitrary polygonal domains using either Euclidean or geodesic distance.
This vignette demonstrates the use of the package on a representative
example with a complex domain boundary.

### Example data

We first define a concave, irregular spatial domain and a set of
generator points with heterogeneous weights.

``` r
library(sf)
library(weightedVoronoi)

crs_use <- 32636  # projected CRS (metres)

boundary_sf <- make_irregular_boundary()

set.seed(42)

pts <- st_sample(boundary_sf, size = 8, type = "random")
points_sf <- st_sf(
  village = paste0("V", seq_along(pts)),
  population = round(exp(rnorm(length(pts), log(300), 0.8))),
  geometry = pts,
  crs = crs_use
)

plot(st_geometry(boundary_sf), border = "black", lwd = 2)
plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")
```

### Weighted Euclidean tessellation

We first generate a weighted tessellation using Euclidean distance.
Generator influence is scaled using a logarithmic transformation of the
population weights.

``` r
out_euc <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = 10,
  weight_transform = log10,
  distance = "euclidean",
  clip_to_boundary = TRUE,
  verbose = FALSE
)

plot(st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_euc$polygons["generator_id"], add = TRUE)
plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")
```

### Weighted geodesic tessellation

We now generate a tessellation using geodesic distance, in which
distances are constrained to paths entirely within the domain.

``` r
out_geo <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = 10,
  weight_transform = log10,
  distance = "geodesic",
  close_mask = TRUE,
  close_iters = 1,
  verbose = FALSE
)

plot(st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_geo$polygons["generator_id"], add = TRUE)
plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")
```

#### Geodesic tessellations with elevation-dependent resistance

In many ecological and social–ecological applications, effective
distance is shaped not only by domain geometry but also by environmental
resistance. Topography, land cover, or infrastructure may increase the
cost of movement between locations. To accommodate this,
`weightedVoronoi` allows geodesic distance calculations to be modified
by an isotropic resistance surface.

When a digital elevation model (DEM) is supplied and
`use_tobler = TRUE`, movement cost is adjusted using Tobler’s hiking
function, which increases effective distance across steep slopes.

``` r
library(sf)
library(terra)
library(weightedVoronoi)

# 1) Use a projected CRS (metres)
crs_use <- 32636

# 2) Define a rectangular boundary (fast + robust for demos)
boundary_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0),
    c(1200, 0),
    c(1200, 800),
    c(0, 800),
    c(0, 0)
  )))),
  crs = crs_use
)

# 3) Sample generator points INSIDE boundary (guaranteed inside)
set.seed(1)
n_pts <- 8
pts <- st_sample(boundary_sf, size = n_pts, type = "random")

points_sf <- st_sf(
  id = seq_len(length(pts)),
  population = round(exp(rnorm(length(pts), log(300), 0.8))),  # skewed weights
  geometry = pts,
  crs = crs_use
)

# Sanity check (should be all TRUE)
inside <- st_within(points_sf, boundary_sf, sparse = FALSE)[, 1]
stopifnot(all(inside))

# 4) Choose resolution (cellsize) — keep consistent everywhere
cellsize <- 10
bnd_v <- terra::vect(boundary_sf)

# 5) Create a synthetic DEM with a strong ridge (easy to see)
dem <- terra::rast(
  ext = terra::ext(bnd_v),
  resolution = cellsize,
  crs = terra::crs(bnd_v)
)

xy <- terra::crds(dem, df = TRUE)

# Ridge parameters (tweak height to strengthen/weaken effect)
y0 <- (min(xy$y) + max(xy$y)) / 2
sigma <- (max(xy$y) - min(xy$y)) / 14
height <- 1000

terra::values(dem) <- height * exp(-((xy$y - y0)^2) / (2 * sigma^2))

par(mfrow = c(1,1))
# Optional plot check
plot(st_geometry(boundary_sf), border = "black", lwd = 2, main = "Demo inputs")
plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")


out_geo_plain <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  verbose = FALSE
)

out_geo_dem <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  dem_rast = dem,
  use_tobler = TRUE,
  verbose = FALSE
)


plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_geo_plain$polygons["generator_id"], main = "Geodesic", add = TRUE)
plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")

plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_geo_dem$polygons["generator_id"], main = "Geodesic + Tobler ridge resistance", add = TRUE)
plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")
```

## Comparing Euclidean and geodesic tessellations

Although both tessellations fully partition the domain, they differ in
how influence propagates through concave boundaries and narrow
corridors. The Euclidean tessellation allocates regions based on
straight-line distance within the rasterised domain, whereas the
geodesic tessellation respects the domain geometry and prevents
influence across inaccessible areas.

These differences are reflected both visually and in the generator-level
summary tables returned by the software.

``` r
out_euc$summary
out_geo$summary
```

## Practical considerations

- Resolution: Smaller raster cell sizes improve boundary fidelity but
  increase computation time.

- Weights: Alternative weight transformations (e.g. square-root or
  power-law) can be supplied to reflect different assumptions about
  generator influence.

- Distance metric: Geodesic distance is recommended when domain geometry
  strongly constrains access or interaction.

## Summary

This vignette illustrates how weightedVoronoi can be used to generate
reproducible, weighted spatial tessellations that respect complex
boundaries and heterogeneous generator influence. The approach provides
a flexible alternative to standard Voronoi methods in ecological and
social–ecological applications.
