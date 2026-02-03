# helper-data.R
# Shared synthetic objects for tests

library(sf)
library(terra)

make_rect_domain <- function(crs_use = 32636) {
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
  boundary_sf
}

make_points_in_domain <- function(boundary_sf, n = 6, seed = 1) {
  set.seed(seed)
  pts <- st_sample(boundary_sf, size = n, type = "random")
  points_sf <- st_sf(
    id = seq_len(length(pts)),
    population = round(exp(rnorm(length(pts), log(300), 0.7))),
    geometry = pts,
    crs = st_crs(boundary_sf)
  )
  # sanity: all inside
  inside <- st_within(points_sf, boundary_sf, sparse = FALSE)[, 1]
  stopifnot(all(inside))
  points_sf
}

make_domain_mask <- function(boundary_sf, res) {
  bnd_v <- terra::vect(boundary_sf)
  r <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
  r <- terra::setValues(r, 1)
  terra::mask(r, bnd_v)
}
