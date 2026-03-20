test_that("temporal tessellation returns allocation stack", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  pts_t1 <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  pts_t2 <- st_sf(
    population = c(2, 1),
    geometry = st_sfc(
      st_point(c(250, 500)),
      st_point(c(750, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_time(
    points_list = list(t1 = pts_t1, t2 = pts_t2),
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    verbose = FALSE
  )
  
  expect_true(inherits(out$allocations, "SpatRaster"))
  expect_equal(terra::nlyr(out$allocations), 2)
  expect_equal(names(out$allocations), c("t1", "t2"))
})

test_that("temporal tessellation can keep summaries", {
  library(sf)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  pts_t1 <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  pts_t2 <- st_sf(
    population = c(2, 1),
    geometry = st_sfc(
      st_point(c(250, 500)),
      st_point(c(750, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_time(
    points_list = list(t1 = pts_t1, t2 = pts_t2),
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    keep_summaries = TRUE,
    verbose = FALSE
  )
  
  expect_true(is.list(out$summaries))
  expect_equal(length(out$summaries), 2)
  expect_true(all(vapply(out$summaries, is.data.frame, logical(1))))
})

test_that("temporal tessellation accepts matching resistance list", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  pts_t1 <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  pts_t2 <- st_sf(
    population = c(2, 1),
    geometry = st_sfc(
      st_point(c(250, 500)),
      st_point(c(750, 500))
    ),
    crs = crs_use
  )
  
  r1 <- terra::rast(ext = terra::ext(0, 1000, 0, 1000), resolution = 50, crs = crs_use)
  terra::values(r1) <- 1
  
  r2 <- terra::rast(ext = terra::ext(0, 1000, 0, 1000), resolution = 50, crs = crs_use)
  terra::values(r2) <- 2
  
  out <- weighted_voronoi_time(
    points_list = list(t1 = pts_t1, t2 = pts_t2),
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    geodesic_engine = "classic",
    resistance_list = list(r1, r2),
    keep_summaries = FALSE,
    weight_model = "additive",
    verbose = FALSE
  )
  
  expect_true(inherits(out$allocations, "SpatRaster"))
  expect_equal(terra::nlyr(out$allocations), 2)
})

test_that("temporal tessellation validates list lengths", {
  library(sf)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  pts_t1 <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  pts_t2 <- st_sf(
    population = c(2, 1),
    geometry = st_sfc(
      st_point(c(250, 500)),
      st_point(c(750, 500))
    ),
    crs = crs_use
  )
  
  expect_error(
    weighted_voronoi_time(
      points_list = list(t1 = pts_t1, t2 = pts_t2),
      weight_col = "population",
      boundary_sf = boundary_sf,
      time_index = "t1",
      verbose = FALSE
    ),
    "time_index must have the same length"
  )
})

test_that("temporal tessellation returns change map and persistence", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:32636"
  
  boundary_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))),
    crs = crs_use
  )
  
  pts_t1 <- st_sf(
    population = c(1, 2),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(800, 500))
    ),
    crs = crs_use
  )
  
  pts_t2 <- st_sf(
    population = c(2, 1),
    geometry = st_sfc(
      st_point(c(250, 500)),
      st_point(c(750, 500))
    ),
    crs = crs_use
  )
  
  out <- weighted_voronoi_time(
    points_list = list(time1 = pts_t1, time2 = pts_t2),
    weight_col = "population",
    boundary_sf = boundary_sf,
    distance = "geodesic",
    geodesic_engine = "multisource",
    weight_model = "additive",
    verbose = FALSE
  )
  
  expect_true(inherits(out$change_map_first_last, "SpatRaster"))
  expect_true(inherits(out$persistence, "SpatRaster"))
  
  change_vals <- terra::values(out$change_map_first_last)
  persist_vals <- terra::values(out$persistence)
  
  change_vals <- change_vals[is.finite(change_vals)]
  persist_vals <- persist_vals[is.finite(persist_vals)]
  
  expect_true(all(change_vals %in% c(0, 1)))
  expect_true(all(persist_vals %in% c(0, 1)))
})