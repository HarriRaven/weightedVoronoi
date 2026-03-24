library(sf)
library(terra)
library(weightedVoronoi)

make_boundary <- function(crs_use = "EPSG:32636", size = 2000) {
  st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(0, 0),
      c(size, 0),
      c(size, size),
      c(0, size),
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

time_expr <- function(expr) {
  t <- system.time(force(expr))
  data.frame(
    elapsed = unname(t["elapsed"]),
    user = unname(t["user.self"]),
    system = unname(t["sys.self"]),
    stringsAsFactors = FALSE
  )
}

benchmark_oneoff <- function(boundary_sf, n_points, res, seed = 1) {
  points_sf <- make_points(boundary_sf, n_points, seed = seed)
  
  out <- list()
  
  out[[1]] <- cbind(
    workflow = "oneoff",
    engine = "classic",
    n_points = n_points,
    res = res,
    n_times = NA,
    n_sim = NA,
    time_expr(
      weighted_voronoi_domain(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        res = res,
        distance = "geodesic",
        geodesic_engine = "classic",
        weight_model = "additive",
        verbose = FALSE
      )
    )
  )
  
  out[[2]] <- cbind(
    workflow = "oneoff",
    engine = "multisource",
    n_points = n_points,
    res = res,
    n_times = NA,
    n_sim = NA,
    time_expr(
      weighted_voronoi_domain(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        res = res,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        verbose = FALSE
      )
    )
  )
  
  ctx <- prepare_geodesic_context(
    boundary_sf = boundary_sf,
    res = res,
    anisotropy = "none",
    geodesic_engine = "multisource"
  )
  
  out[[3]] <- cbind(
    workflow = "oneoff",
    engine = "prepared_multisource",
    n_points = n_points,
    res = res,
    n_times = NA,
    n_sim = NA,
    time_expr(
      weighted_voronoi_domain(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        res = res,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        prepared = ctx,
        verbose = FALSE
      )
    )
  )
  
  do.call(rbind, out)
}

benchmark_temporal <- function(boundary_sf, n_points, res, n_times, seed = 1) {
  points_list <- lapply(seq_len(n_times), function(i) {
    make_points(boundary_sf, n_points, seed = seed + i)
  })
  names(points_list) <- paste0("t", seq_len(n_times))
  
  cbind(
    workflow = "temporal",
    engine = "multisource",
    n_points = n_points,
    res = res,
    n_times = n_times,
    n_sim = NA,
    time_expr(
      weighted_voronoi_time(
        points_list = points_list,
        weight_col = "population",
        boundary_sf = boundary_sf,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        res = res,
        verbose = FALSE
      )
    )
  )
}

benchmark_uncertainty <- function(boundary_sf, n_points, res, n_sim, seed = 1) {
  points_sf <- make_points(boundary_sf, n_points, seed = seed)
  
  cbind(
    workflow = "uncertainty",
    engine = "multisource",
    n_points = n_points,
    res = res,
    n_times = NA,
    n_sim = n_sim,
    time_expr(
      weighted_voronoi_uncertainty(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        n_sim = n_sim,
        weight_sd = 0.3,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        res = res,
        verbose = FALSE,
        warn_zero_entropy = FALSE,
        seed = seed
      )
    )
  )
}

boundary_sf <- make_boundary()

results <- list()

# One-off grid
oneoff_grid <- expand.grid(
  n_points = c(5, 10, 20, 40, 80),
  res = c(100, 50, 25),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(oneoff_grid))) {
  g <- oneoff_grid[i, ]
  message(sprintf("One-off: n_points=%d, res=%d", g$n_points, g$res))
  results[[length(results) + 1]] <- benchmark_oneoff(
    boundary_sf = boundary_sf,
    n_points = g$n_points,
    res = g$res,
    seed = 100 + i
  )
}

# Temporal grid
temporal_grid <- expand.grid(
  n_points = c(10, 20, 40),
  res = c(100, 50),
  n_times = c(2, 4, 8),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(temporal_grid))) {
  g <- temporal_grid[i, ]
  message(sprintf("Temporal: n_points=%d, res=%d, n_times=%d", g$n_points, g$res, g$n_times))
  results[[length(results) + 1]] <- benchmark_temporal(
    boundary_sf = boundary_sf,
    n_points = g$n_points,
    res = g$res,
    n_times = g$n_times,
    seed = 1000 + i
  )
}

# Uncertainty grid
uncertainty_grid <- expand.grid(
  n_points = c(10, 20, 40),
  res = c(100, 50),
  n_sim = c(10, 25, 50, 100),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(uncertainty_grid))) {
  g <- uncertainty_grid[i, ]
  message(sprintf("Uncertainty: n_points=%d, res=%d, n_sim=%d", g$n_points, g$res, g$n_sim))
  results[[length(results) + 1]] <- benchmark_uncertainty(
    boundary_sf = boundary_sf,
    n_points = g$n_points,
    res = g$res,
    n_sim = g$n_sim,
    seed = 2000 + i
  )
}

bench_df <- do.call(rbind, results)
rownames(bench_df) <- NULL

dir.create("bench", showWarnings = FALSE, recursive = TRUE)
write.csv(bench_df, "bench/benchmark_results.csv", row.names = FALSE)

print(bench_df)