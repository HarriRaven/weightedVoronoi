test_that("Same inputs produce identical allocation (deterministic)", {
  boundary_sf <- make_rect_domain()
  points_sf <- make_points_in_domain(boundary_sf, n = 6, seed = 3)
  res <- 25
  
  out1 <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "euclidean",
    verbose = FALSE
  )
  
  out2 <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "euclidean",
    verbose = FALSE
  )
  
  v1 <- terra::values(out1$allocation, mat = FALSE)
  v2 <- terra::values(out2$allocation, mat = FALSE)
  expect_identical(v1, v2)
})
