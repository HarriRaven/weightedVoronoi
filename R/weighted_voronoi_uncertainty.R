#' Uncertainty-aware weighted tessellation
#'
#' Repeats weighted tessellation under stochastic perturbation of generator
#' weights and summarises the results as per-cell membership probabilities,
#' modal allocation, and entropy.
#'
#' @param points_sf An `sf` POINT object containing generator locations and attributes.
#' @param weight_col Character. Name of the weight column in `points_sf`.
#' @param boundary_sf An `sf` POLYGON/MULTIPOLYGON defining the domain.
#' @param n_sim Integer. Number of simulation runs.
#' @param weight_sd Optional numeric. Standard deviation of lognormal weight
#'   perturbation on the log scale. If `NULL`, no perturbation is applied and the
#'   same tessellation is repeated.
#' @param distance Character. One of `"euclidean"` or `"geodesic"`.
#' @param geodesic_engine Character. Geodesic engine to use when
#'   `distance = "geodesic"`. Defaults to `"multisource"`.
#' @param keep_simulations Logical. If `TRUE`, return the simulated allocation
#'   rasters as a stack.
#' @param seed Optional integer random seed for reproducibility.
#' @param verbose Logical. If `TRUE`, prints progress.
#' @param ... Additional arguments passed to [weighted_voronoi_domain()].
#' @param warn_zero_entropy Logical. If `TRUE`, warn when all entropy values are
#'   zero across the domain.
#' @param res Numeric. Raster resolution in CRS units (e.g. metres).
#' @param prepared Optional prepared geodesic context created by
#'   [prepare_geodesic_context()] for repeated compatible geodesic runs.
#'
#' @details
#' This first implementation supports uncertainty in generator weights only.
#' Weights are perturbed independently across simulations using a lognormal
#' multiplicative model:
#'
#' `w_sim = w * exp(N(0, weight_sd))`
#'
#' The output includes:
#' \describe{
#'   \item{probabilities}{A `terra::SpatRaster` with one layer per generator,
#'   containing the fraction of simulations in which each cell was assigned to
#'   that generator.}
#'   \item{modal_allocation}{A `terra::SpatRaster` giving the most probable
#'   generator for each cell.}
#'   \item{entropy}{A `terra::SpatRaster` showing per-cell uncertainty. Higher
#'   values indicate less stable allocation across simulations.}
#' }
#'
#' @return A list with probability surfaces, modal allocation, entropy, and
#'   optionally the full simulation stack.
#' @export
weighted_voronoi_uncertainty <- function(
    points_sf,
    weight_col,
    boundary_sf,
    n_sim = 100,
    weight_sd = NULL,
    distance = c("euclidean", "geodesic"),
    geodesic_engine = c("multisource", "classic"),
    res = 20,
    keep_simulations = FALSE,
    seed = NULL,
    warn_zero_entropy = TRUE,
    prepared = NULL,
    verbose = TRUE,
    ...
) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  
  distance <- match.arg(distance)
  geodesic_engine <- match.arg(geodesic_engine)
  
  
  if (!inherits(points_sf, "sf")) stop("points_sf must be an sf object.")
  if (!inherits(boundary_sf, "sf")) stop("boundary_sf must be an sf object.")
  if (!weight_col %in% names(points_sf)) stop("weight_col not found in points_sf.")
  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 1) {
    stop("n_sim must be a finite integer >= 1.")
  }
  n_sim <- as.integer(n_sim)
  
  if (!is.null(weight_sd) && (!is.numeric(weight_sd) || length(weight_sd) != 1 ||
                              !is.finite(weight_sd) || weight_sd < 0)) {
    stop("weight_sd must be NULL or a finite numeric >= 0.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  dots <- list(...)
  dot1 <- function(name, default) {
    if (name %in% names(dots)) dots[[name]] else default
  }
  
  base_weights <- as.numeric(points_sf[[weight_col]])
  if (any(!is.finite(base_weights)) || any(base_weights <= 0)) {
    stop("Input weights must be finite and > 0 before uncertainty perturbation.")
  }
  
  gen_ids <- seq_len(nrow(points_sf))
  n_gen <- length(gen_ids)
  
  prepared_ctx <- NULL
  if (distance == "geodesic") {
    if (!is.null(prepared)) {
      prepared_ctx <- prepared
    } else {
      prepared_ctx <- prepare_geodesic_context(
        boundary_sf = boundary_sf,
        res = res,
        close_mask = dot1("close_mask", TRUE),
        close_iters = dot1("close_iters", 1),
        resistance_rast = dot1("resistance_rast", NULL),
        dem_rast = dot1("dem_rast", NULL),
        use_tobler = dot1("use_tobler", TRUE),
        tobler_v0_kmh = dot1("tobler_v0_kmh", 6),
        tobler_a = dot1("tobler_a", 3.5),
        tobler_b = dot1("tobler_b", 0.05),
        min_speed_kmh = dot1("min_speed_kmh", 0.25),
        anisotropy = dot1("anisotropy", "none"),
        uphill_factor = dot1("uphill_factor", 1),
        downhill_factor = dot1("downhill_factor", 1),
        geodesic_engine = geodesic_engine
      )
    }
  }
  
  allocation_list <- if (keep_simulations) vector("list", n_sim) else NULL
  
  template_rast <- NULL
  inside_domain <- NULL
  count_mat <- NULL  # ncell x n_gen integer matrix
  
  for (i in seq_len(n_sim)) {
    if (verbose && (i == 1 || i == n_sim || i %% 10 == 0)) {
      message(sprintf("Uncertainty simulation %d / %d", i, n_sim))
    }
    
    points_i <- points_sf
    
    if (!is.null(weight_sd) && weight_sd > 0) {
      w_sim <- base_weights * exp(stats::rnorm(length(base_weights), mean = 0, sd = weight_sd))
      points_i[[weight_col]] <- w_sim
    }
    
    if (distance == "geodesic") {
      out_i <- weighted_voronoi_geodesic(
        points_sf = points_i,
        weight_col = weight_col,
        boundary_sf = boundary_sf,
        geodesic_engine = geodesic_engine,
        res = res,
        prepared = prepared_ctx,
        return_polygons = FALSE,
        verbose = FALSE,
        ...
      )
      allocation_i <- out_i$allocation
      
    } else {
      out_i <- weighted_voronoi_domain(
        points_sf = points_i,
        weight_col = weight_col,
        boundary_sf = boundary_sf,
        distance = distance,
        geodesic_engine = geodesic_engine,
        res = res,
        verbose = FALSE,
        ...
      )
      allocation_i <- out_i$allocation
    }
    
    alloc_vals <- terra::values(allocation_i, mat = FALSE)
    
    if (is.null(template_rast)) {
      template_rast <- allocation_i
      inside_domain <- !is.na(alloc_vals)
      n_cells <- length(alloc_vals)
      count_mat <- matrix(0L, nrow = n_cells, ncol = n_gen)
    }
    
    for (j in seq_len(n_gen)) {
      count_mat[, j] <- count_mat[, j] + as.integer(!is.na(alloc_vals) & alloc_vals == gen_ids[j])
    }
    
    if (keep_simulations) {
      names(allocation_i) <- paste0("sim_", i)
      allocation_list[[i]] <- allocation_i
    }
  }
  
  # Build probability rasters from count matrix
  prob_list <- vector("list", n_gen)
  for (j in seq_len(n_gen)) {
    vals <- rep(NA_real_, nrow(count_mat))
    vals[inside_domain] <- count_mat[inside_domain, j] / n_sim
    
    r <- template_rast
    r <- terra::setValues(r, vals)
    names(r) <- paste0("gen_", gen_ids[j])
    prob_list[[j]] <- r
  }
  
  prob_stack <- terra::rast(prob_list)
  
  modal_allocation <- terra::app(prob_stack, fun = function(x) {
    if (all(is.na(x))) return(NA_real_)
    which.max(x)
  })
  names(modal_allocation) <- "modal_allocation"
  
  entropy <- terra::app(prob_stack, fun = function(p) {
    p <- p[is.finite(p) & p > 0]
    if (!length(p)) return(NA_real_)
    -sum(p * log(p))
  })
  names(entropy) <- "entropy"
  
  ent_vals <- terra::values(entropy, mat = FALSE)
  ent_vals <- ent_vals[is.finite(ent_vals)]
  
  if (warn_zero_entropy && length(ent_vals) && all(ent_vals == 0)) {
    warning(
      "All entropy values are zero. Simulated perturbations did not change cell allocation. ",
      "Consider increasing weight_sd or using a setup where weight effects are larger relative to spatial distances."
    )
  }
  
  out <- list(
    probabilities = prob_stack,
    modal_allocation = modal_allocation,
    entropy = entropy,
    n_sim = n_sim,
    weight_sd = weight_sd,
    distance = distance,
    geodesic_engine = geodesic_engine
  )
  
  if (keep_simulations) {
    alloc_stack <- terra::rast(allocation_list)
    out$simulations <- alloc_stack
  }
  
  out
}