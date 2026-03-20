library(sf)
library(terra)
library(weightedVoronoi)

make_boundary <- function(crs_use = "EPSG:32636") {
  st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0, 0),
      c(2000, 0),
      c(2000, 2000),
      c(0, 2000),
      c(0, 0)
    )))),
    crs = crs_use
  )
}

make_points <- function(boundary_sf, n, seed = 1) {
  set.seed(seed)
  pts <- st_sample(boundary_sf, size = n, type = "random")
  st_sf(
    population = sample(1:5, length(pts), replace = TRUE),
    geometry = pts,
    crs = st_crs(boundary_sf)
  )
}

run_case <- function(n_points, res, engine) {
  boundary_sf <- make_boundary()
  points_sf <- make_points(boundary_sf, n_points, seed = 100 + n_points)
  
  t <- system.time({
    out <- weighted_voronoi_domain(
      points_sf = points_sf,
      weight_col = "population",
      boundary_sf = boundary_sf,
      res = res,
      distance = "geodesic",
      weight_model = "additive",
      geodesic_engine = engine,
      verbose = FALSE
    )
  })
  
  data.frame(
    n_points = n_points,
    res = res,
    engine = engine,
    elapsed = unname(t["elapsed"]),
    user = unname(t["user.self"]),
    system = unname(t["sys.self"]),
    stringsAsFactors = FALSE
  )
}

grid <- expand.grid(
  n_points = c(5, 10, 20, 40),
  res = c(100, 50),
  engine = c("classic", "multisource"),
  stringsAsFactors = FALSE
)

results <- do.call(
  rbind,
  lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    message(sprintf("Running n_points=%d, res=%d, engine=%s",
                    g$n_points, g$res, g$engine))
    run_case(g$n_points, g$res, g$engine)
  })
)

print(results)