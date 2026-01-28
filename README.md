
<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

# weightedVoronoi

Tools for weighted spatial tessellation using Euclidean and geodesic
distances within constrained polygons. Produces complete, connected
partitions that respect complex boundaries and heterogeneous point
weights.

## Installation

\`\`\`r \# install.packages(“remotes”)
remotes::install_github(“HarriRaven/weightedVoronoi”)

library(sf) library(weightedVoronoi)

# Use a projected CRS (units in metres)

crs_use \<- 32636

# Domain polygon (simple rectangle for speed)

boundary_sf \<- st_sf( geometry = st_sfc(st_polygon(list(rbind( c(0, 0),
c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0) )))), crs = crs_use )

# Generator points with weights

points_sf \<- st_sf( village = paste0(“V”, 1:5), population = c(50, 200,
1000, 150, 400), geometry = st_sfc( st_point(c(200, 200)),
st_point(c(800, 250)), st_point(c(500, 500)), st_point(c(250, 800)),
st_point(c(750, 750)) ), crs = crs_use )

# Weighted Euclidean tessellation

out_euc \<- weighted_voronoi_domain( points_sf = points_sf, weight_col =
“population”, boundary_sf = boundary_sf, res = 20, weight_transform =
log10, distance = “euclidean”, verbose = FALSE )

# Weighted geodesic tessellation (domain-constrained shortest path distance)

out_geo \<- weighted_voronoi_domain( points_sf = points_sf, weight_col =
“population”, boundary_sf = boundary_sf, res = 20, weight_transform =
log10, distance = “geodesic”, close_mask = TRUE, close_iters = 1,
verbose = FALSE )

# Inspect outputs

names(out_euc) head(out_euc$summary)
out_euc$diagnostics

Outputs

weighted_voronoi_domain() returns:

- polygons: sf object with one polygon per generator (and attributes)

- allocation: terra::SpatRaster assigning each raster cell to a
  generator

- summary: generator-level summary table (area, share, weights, etc.)

- diagnostics: diagnostics and settings (coverage, unreachable fraction
  for geodesic, etc.)

Notes

- Inputs must be in a projected CRS with metric units (e.g. metres).

- res controls the raster resolution and therefore the trade-off between
  speed and boundary fidelity.

- Geodesic tessellations are typically slower than Euclidean
  tessellations because shortest-path distances are computed within the
  domain.

Citation

If you use weightedVoronoi, please cite the associated software note.
