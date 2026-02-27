library(sf)

make_square <- function(size = 1000, crs = 32630) {
  st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(size,0), c(size,size), c(0,size), c(0,0)
    )))),
    crs = crs
  )
}

make_points <- function(crs = 32630) {
  st_sf(
    w = c(1, 2, 3),
    geometry = st_sfc(
      st_point(c(200, 200)),
      st_point(c(800, 250)),
      st_point(c(500, 800))
    ),
    crs = crs
  )
}

#-------------------------------------------------------#
#Test A - weight model equivalence----
#-------------------------------------------------------#

test_that("weight_model power with weight_power=1 matches multiplicative", {
  b <- make_square()
  p <- make_points()
  
  out1 <- weighted_voronoi_domain(
    p, "w", b, res = 50,
    distance = "euclidean",
    weight_model = "multiplicative",
    verbose = FALSE
  )
  
  out2 <- weighted_voronoi_domain(
    p, "w", b, res = 50,
    distance = "euclidean",
    weight_model = "power",
    weight_power = 1,
    verbose = FALSE
  )
  
  v1 <- terra::values(out1$allocation)
  v2 <- terra::values(out2$allocation)
  
  expect_identical(v1, v2)
})

#-------------------------------------------------------#
#Test B - compose_resistance() aligns and returns >0----
#-------------------------------------------------------#

test_that("compose_resistance aligns to template and stays >0", {
  b <- make_square()
  p <- make_points()
  
  # Build a template grid by running a tiny euclidean tessellation and grabbing allocation grid
  tmp <- weighted_voronoi_domain(p, "w", b, res = 50, distance="euclidean", verbose=FALSE)$allocation
  
  r1 <- tmp; terra::values(r1) <- 2
  r2 <- tmp; terra::values(r2) <- 3
  
  R <- compose_resistance(r1, r2, template = tmp, method = "multiply")
  expect_true(inherits(R, "SpatRaster"))
  
  # Check strictly > 0 where not NA
  vals <- terra::values(R)
  vals <- vals[is.finite(vals)]
  expect_true(all(vals > 0))
})

#-------------------------------------------------------#
#Test C - geodesic uses resistance_rast (precedence)----
#-------------------------------------------------------#

test_that("geodesic reacts to resistance_rast (precedence works)", {
  b <- make_square()
  p <- make_points()
  
  # Baseline: uniform resistance (default)
  out_base <- weighted_voronoi_domain(
    p, "w", b, res = 50,
    distance = "geodesic",
    close_mask = FALSE,
    verbose = FALSE
  )
  
  # Now make a high-friction vertical band in the middle
  # Build a template raster with same grid by using an euclidean allocation grid
  tmp <- weighted_voronoi_domain(p, "w", b, res = 50, distance="euclidean", verbose=FALSE)$allocation
  R <- tmp; terra::values(R) <- 1
  
  xy <- terra::xyFromCell(R, 1:terra::ncell(R))
  mid_band <- xy[,1] > 450 & xy[,1] < 550
  vals <- terra::values(R)
  vals[mid_band] <- 50
  terra::values(R) <- vals
  
  out_res <- weighted_voronoi_domain(
    p, "w", b, res = 50,
    distance = "geodesic",
    resistance_rast = R,
    close_mask = FALSE,
    verbose = FALSE
  )
  
  # Allocation should differ in at least some cells
  v1 <- terra::values(out_base$allocation)
  v2 <- terra::values(out_res$allocation)
  expect_true(any(v1 != v2, na.rm = TRUE))
})

#-------------------------------------------------------#
#Test D - impermeable barrier doesn’t crash and changes allocation----
#-------------------------------------------------------#

test_that("impermeable barrier works in geodesic via user-side preprocessing", {
  b <- make_square()
  p <- make_points()
  
  tmp <- weighted_voronoi_domain(p, "w", b, res = 50, distance="euclidean", verbose=FALSE)$allocation
  R <- tmp; terra::values(R) <- 1
  
  # Create a vertical line barrier x=500
  river <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_linestring(rbind(c(500,0), c(500,1000)))),
    crs = sf::st_crs(b)
  )
  
  Rb <- add_barriers(R, river, permeability="impermeable", width=25)
  
  out0 <- weighted_voronoi_domain(p, "w", b, res=50, distance="geodesic",
                                  resistance_rast=R, close_mask=FALSE, verbose=FALSE)
  out1 <- weighted_voronoi_domain(p, "w", b, res=50, distance="geodesic",
                                  resistance_rast=Rb, close_mask=FALSE, verbose=FALSE)
  
  v0 <- terra::values(out0$allocation)
  v1 <- terra::values(out1$allocation)
  
  expect_true(any(v0 != v1, na.rm = TRUE))
})