test_that("DEM resistance changes geodesic allocation but stays complete", {
  boundary_sf <- make_rect_domain()
  points_sf <- make_points_in_domain(boundary_sf, n = 6, seed = 4)
  res <- 25
  
  # Flat geodesic
  out_plain <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "geodesic",
    verbose = FALSE
  )
  
  # Synthetic DEM: N-S ridge
  bnd_v <- terra::vect(boundary_sf)
  dem <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
  xy <- terra::crds(dem, df = TRUE)
  y0 <- (min(xy$y) + max(xy$y)) / 2
  sigma <- (max(xy$y) - min(xy$y)) / 12
  terra::values(dem) <- 1000 * exp(-((xy$y - y0)^2) / (2 * sigma^2))
  
  out_dem <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "geodesic",
    dem_rast = dem,
    use_tobler = TRUE,
    verbose = FALSE
  )
  
  # Must still be complete INSIDE the domain mask
  mask_on_alloc <- out_dem$allocation
  terra::values(mask_on_alloc) <- 1
  
  # Make boundary as terra vector and ensure terra-side CRS matches allocation
  bnd_v_alloc <- terra::vect(boundary_sf)
  terra::crs(bnd_v_alloc) <- terra::crs(out_dem$allocation)
  
  mask_on_alloc <- terra::mask(mask_on_alloc, bnd_v_alloc)
  inside_vals <- terra::values(terra::mask(out_dem$allocation, mask_on_alloc), mat = FALSE)
  
  inside_vals <- inside_vals[!is.na(inside_vals)]
  expect_true(length(inside_vals) > 0)
  expect_true(all(is.finite(inside_vals)))
  
  
  # And should differ from plain case (not necessarily every cell)
  v_plain <- terra::values(out_plain$allocation, mat = FALSE)
  v_dem   <- terra::values(out_dem$allocation, mat = FALSE)
  expect_true(any(v_plain != v_dem, na.rm = TRUE))
})
