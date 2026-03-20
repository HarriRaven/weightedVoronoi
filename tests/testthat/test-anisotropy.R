test_that("terrain anisotropy runs and returns complete allocation", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:3857"
  
  boundary_sf <- st_sf(
    id = 1,
    geometry = st_sfc(
      st_polygon(list(rbind(
        c(0, 0),
        c(1000, 0),
        c(1000, 1000),
        c(0, 1000),
        c(0, 0)
      ))),
      crs = crs_use
    )
  )
  
  points_sf <- st_sf(
    population = c(1, 1),
    geometry = st_sfc(
      st_point(c(200, 500)),
      st_point(c(750, 500)),
      crs = crs_use
    )
  )
  
  dem_rast <- rast(
    ext = ext(0, 1000, 0, 1000),
    resolution = 50,
    crs = crs_use
  )
  
  xy <- crds(dem_rast, df = TRUE)
  values(dem_rast) <- xy$x * 10
  
  out <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = 50,
    distance = "geodesic",
    dem_rast = dem_rast,
    anisotropy = "terrain",
    uphill_factor = 3,
    downhill_factor = 1.25,
    verbose = FALSE
  )
  
  expect_true(inherits(out$allocation, "SpatRaster"))
  expect_true(inherits(out$polygons, "sf"))
  expect_equal(nrow(out$summary), 2)
  
  vals <- values(out$allocation)
  expect_true(any(!is.na(vals)))
  expect_false(all(is.na(vals)))
})

test_that("terrain anisotropy is accepted and reported in diagnostics", {
  library(sf)
  library(terra)
  
  crs_use <- "EPSG:3857"
  
  boundary_sf <- st_sf(
    id = 1,
    geometry = st_sfc(
      st_polygon(list(rbind(
        c(0, 0),
        c(1000, 0),
        c(1000, 1000),
        c(0, 1000),
        c(0, 0)
      ))),
      crs = crs_use
    )
  )
  
  points_sf <- st_sf(
    population = c(1, 1),
    geometry = st_sfc(
      st_point(c(150, 500)),
      st_point(c(800, 500)),
      crs = crs_use
    )
  )
  
  dem_rast <- rast(
    ext = ext(0, 1000, 0, 1000),
    resolution = 25,
    crs = crs_use
  )
  
  xy <- crds(dem_rast, df = TRUE)
  values(dem_rast) <- xy$x * 50
  
  out_aniso <- weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = 25,
    distance = "geodesic",
    dem_rast = dem_rast,
    anisotropy = "terrain",
    uphill_factor = 6,
    downhill_factor = 1.5,
    verbose = FALSE
  )
  
  expect_equal(out_aniso$diagnostics$anisotropy, "terrain")
  expect_equal(out_aniso$diagnostics$uphill_factor, 6)
  expect_equal(out_aniso$diagnostics$downhill_factor, 1.5)
  expect_true(inherits(out_aniso$allocation, "SpatRaster"))
  expect_true(inherits(out_aniso$polygons, "sf"))
})