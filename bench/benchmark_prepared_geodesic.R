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

time_expr <- function(expr) {
  t <- system.time(force(expr))
  data.frame(
    elapsed = unname(t["elapsed"]),
    user = unname(t["user.self"]),
    system = unname(t["sys.self"]),
    stringsAsFactors = FALSE
  )
}

boundary_sf <- make_boundary()
points_sf <- make_points(boundary_sf, n = 40, seed = 123)

cat("=== One-off geodesic benchmark ===\n")

res_oneoff <- rbind(
  cbind(
    scenario = "classic",
    time_expr(
      weighted_voronoi_domain(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        res = 40,
        distance = "geodesic",
        geodesic_engine = "classic",
        weight_model = "additive",
        verbose = FALSE
      )
    )
  ),
  cbind(
    scenario = "multisource",
    time_expr(
      weighted_voronoi_domain(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        res = 40,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        verbose = FALSE
      )
    )
  ),
  cbind(
    scenario = "prepared_multisource",
    {
      ctx <- prepare_geodesic_context(
        boundary_sf = boundary_sf,
        res = 40,
        anisotropy = "none",
        geodesic_engine = "multisource"
      )
      time_expr(
        weighted_voronoi_domain(
          points_sf = points_sf,
          weight_col = "population",
          boundary_sf = boundary_sf,
          res = 40,
          distance = "geodesic",
          geodesic_engine = "multisource",
          weight_model = "additive",
          prepared = ctx,
          verbose = FALSE
        )
      )
    }
  )
)

print(res_oneoff)

cat("\n=== Uncertainty benchmark ===\n")

res_uncertainty <- rbind(
  cbind(
    scenario = "uncertainty_multisource",
    time_expr(
      weighted_voronoi_uncertainty(
        points_sf = points_sf,
        weight_col = "population",
        boundary_sf = boundary_sf,
        n_sim = 50,
        weight_sd = 0.3,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        verbose = FALSE,
        warn_zero_entropy = FALSE,
        seed = 1,
        res = 40
      )
    )
  )
)

print(res_uncertainty)

cat("\n=== Temporal benchmark ===\n")

points_list <- list(
  t1 = make_points(boundary_sf, n = 40, seed = 201),
  t2 = make_points(boundary_sf, n = 40, seed = 202),
  t3 = make_points(boundary_sf, n = 40, seed = 203),
  t4 = make_points(boundary_sf, n = 40, seed = 204)
)

res_temporal <- rbind(
  cbind(
    scenario = "temporal_multisource",
    time_expr(
      weighted_voronoi_time(
        points_list = points_list,
        weight_col = "population",
        boundary_sf = boundary_sf,
        distance = "geodesic",
        geodesic_engine = "multisource",
        weight_model = "additive",
        verbose = FALSE,
        res = 40
      )
    )
  )
)

print(res_temporal)