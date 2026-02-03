test_that("Geodesic tessellation runs and yields complete partition", {
  boundary_sf <- make_rect_domain()
  points_sf <- make_points_in_domain(boundary_sf, n = 6, seed = 2)
  res <- 25
  
  out <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "geodesic",
    close_mask = TRUE,
    close_iters = 1,
    verbose = FALSE
  )
  
  expect_s3_class(out$polygons, "sf")
  expect_true(inherits(out$allocation,"SpatRaster"))
  
  # Complete partition: no NA inside mask
  mask <- make_domain_mask(boundary_sf, res = res)
  alloc <- terra::resample(out$allocation, mask, method = "near")
  inside_vals <- terra::values(terra::mask(alloc, mask), mat = FALSE)
  expect_false(any(is.na(inside_vals)))
  
  # If you return unreachable raster/logical, ensure it has correct shape
  if ("unreachable" %in% names(out)) {
    # unreachable could be SpatRaster or logical SpatRaster; just check compatible
    expect_true(inherits(out$unreachable,"SpatRaster"))
  }
})
