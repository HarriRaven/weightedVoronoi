test_that("uncertainty wrapper runs and returns expected structure", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  points_sf <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_uncertainty(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    n_sim = 5,
    weight_sd = 0.1,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    warn_zero_entropy = FALSE,
    verbose = FALSE,
    seed = 123
  )
  
  expect_true(inherits(out$probabilities, "SpatRaster"))
  expect_true(inherits(out$modal_allocation, "SpatRaster"))
  expect_true(inherits(out$entropy, "SpatRaster"))
  expect_equal(out$n_sim, 5)
})

test_that("probabilities sum to one within allocated cells", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  points_sf <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_uncertainty(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    n_sim = 5,
    weight_sd = 0.1,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    warn_zero_entropy = FALSE,
    verbose = FALSE,
    seed = 123
  )
  
  s <- terra::app(out$probabilities, sum)
  vals <- terra::values(s)
  vals <- vals[is.finite(vals)]
  
  expect_true(all(abs(vals - 1) < 1e-6))
})

test_that("deterministic uncertainty run gives crisp probabilities", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  points_sf <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_uncertainty(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    n_sim = 3,
    weight_sd = NULL,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    warn_zero_entropy = FALSE,
    verbose = FALSE,
    seed = 123
  )
  
  vals <- terra::values(out$probabilities)
  vals <- as.vector(vals)
  vals <- vals[is.finite(vals)]
  
  expect_true(all(vals %in% c(0, 1)))
})

test_that("weight uncertainty can generate non-zero entropy in a sensitive setup", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(200,0), c(200,200), c(0,200), c(0,0)
    )))),
    crs = crs_use
  )
  
  points_sf <- st_sf(
    population = c(0.01, 0.02),
    geometry = st_sfc(
      st_point(c(60, 100)),
      st_point(c(140, 100))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_uncertainty(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    n_sim = 50,
    weight_sd = 0.8,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    warn_zero_entropy = FALSE,
    verbose = FALSE,
    seed = 1
  )
  
  ent_vals <- terra::values(out$entropy)
  ent_vals <- ent_vals[is.finite(ent_vals)]
  
  expect_true(any(ent_vals > 0))
})