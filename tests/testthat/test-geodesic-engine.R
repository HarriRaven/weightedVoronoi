test_that("classic geodesic engine still runs", {
  library(sf)
  
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
  
  out <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    geodesic_engine = "classic",
    verbose = FALSE
  )
  
  expect_true(inherits(out$allocation, "SpatRaster"))
})

test_that("multisource geodesic engine runs for additive isotropic geodesics", {
  library(sf)
  
  crs_use <- 32636
  
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
  
  out <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    weight_model = "additive",
    geodesic_engine = "multisource",
    verbose = FALSE
  )
  
  expect_true(inherits(out$allocation, "SpatRaster"))
  expect_true(inherits(out$polygons, "sf"))
})

test_that("multisource engine rejects unsupported settings", {
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
  
  dem_rast <- rast(ext = ext(0, 1000, 0, 1000), resolution = 50, crs = crs_use)
  values(dem_rast) <- 1:ncell(dem_rast)
  
  expect_error(
    weighted_voronoi_domain(
      points_sf = points_sf,
      weight_col = "population",
      boundary_sf = boundary_sf,
      distance = "geodesic",
      weight_model = "multiplicative",
      geodesic_engine = "multisource",
      verbose = FALSE
    ),
    "requires weight_model = 'additive'"
  )
  
  expect_error(
    weighted_voronoi_domain(
      points_sf = points_sf,
      weight_col = "population",
      boundary_sf = boundary_sf,
      distance = "geodesic",
      weight_model = "additive",
      geodesic_engine = "multisource",
      anisotropy = "terrain",
      dem_rast = dem_rast,
      verbose = FALSE
    ),
    "requires anisotropy = 'none'"
  )
})

test_that("multisource matches classic for additive isotropic geodesics", {
  library(sf)
  
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
  
  out_classic <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    weight_model = "additive",
    geodesic_engine = "classic",
    verbose = FALSE
  )
  
  out_multi <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    weight_model = "additive",
    geodesic_engine = "multisource",
    verbose = FALSE
  )
  
  expect_equal(
    unname(as.vector(terra::values(out_classic$allocation))),
    unname(as.vector(terra::values(out_multi$allocation)))
  )
  
  expect_equal(out_classic$summary$n_cells, out_multi$summary$n_cells)
  expect_equal(out_classic$summary$area_m2, out_multi$summary$area_m2)
})