test_that("Euclidean tessellation returns expected structure and complete partition", {
  boundary_sf <- make_rect_domain()
  points_sf <- make_points_in_domain(boundary_sf, n = 6, seed = 1)
  res <- 25
  
  out <- weightedVoronoi::weighted_voronoi_domain(
    points_sf = points_sf,
    weight_col = "population",
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = log10,
    distance = "euclidean",
    clip_to_boundary = TRUE,
    verbose = FALSE
  )
  
  # Structure
  expect_type(out, "list")
  expect_true(all(c("polygons", "allocation", "summary", "diagnostics") %in% names(out)))
  
  expect_s3_class(out$polygons, "sf")
  expect_true(inherits(out$allocation,"SpatRaster"))
  expect_true(is.data.frame(out$summary))
  
  # generator_id present and has correct range
  expect_true("generator_id" %in% names(out$polygons))
  expect_true(all(sort(unique(out$polygons$generator_id)) %in% seq_len(nrow(points_sf))))
  
  # Complete partition: no NA inside domain mask
  mask <- make_domain_mask(boundary_sf, res = res)
  alloc <- out$allocation
  
  # Align in case terra adjusts extents slightly
  alloc <- terra::resample(alloc, mask, method = "near")
  
  inside_vals <- terra::values(terra::mask(alloc, mask), mat = FALSE)
  expect_false(any(is.na(inside_vals)))
  
  # Areas: union polygons ~ boundary area (tolerance due to rasterization)
  b_area <- as.numeric(sf::st_area(sf::st_union(boundary_sf)))
  p_area <- sum(as.numeric(sf::st_area(out$polygons)))
  expect_lt(abs(b_area - p_area) / b_area, 0.05)
})
