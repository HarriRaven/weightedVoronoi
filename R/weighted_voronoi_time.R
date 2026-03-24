#' Temporal weighted tessellation
#'
#' Runs weighted tessellation across a sequence of time-specific point datasets
#' and returns a stack of allocation rasters, with optional polygons and summaries
#' per time step.
#'
#' @param points_list A non-empty list of `sf` POINT objects, one per time step.
#' @param weight_col Character. Name of the weight column present in each element
#'   of `points_list`.
#' @param boundary_sf An `sf` POLYGON/MULTIPOLYGON defining the domain.
#' @param time_index Optional character vector of time labels. Defaults to
#'   `names(points_list)` if present, otherwise `"t1"`, `"t2"`, etc.
#' @param distance Character. One of `"euclidean"` or `"geodesic"`.
#' @param geodesic_engine Character. Geodesic engine to use when
#'   `distance = "geodesic"`.
#' @param resistance_list Optional list of resistance rasters, either length 1
#'   (reused for all times) or the same length as `points_list`.
#' @param dem_list Optional list of DEM rasters, either length 1 (reused for all
#'   times) or the same length as `points_list`.
#' @param keep_polygons Logical. If `TRUE`, return polygons for each time step.
#' @param keep_summaries Logical. If `TRUE`, return summaries for each time step.
#' @param verbose Logical. If `TRUE`, prints progress.
#' @param ... Additional arguments passed to [weighted_voronoi_domain()].
#' @param res Numeric. Raster resolution in CRS units (e.g. metres).
#'
#' @details
#' This first implementation assumes a static boundary and runs each time step
#' independently. Time-varying weights, point locations, and resistance/DEM
#' surfaces are supported by supplying separate inputs for each time step.
#'
#' @return A list containing:
#' \describe{
#'   \item{allocations}{A `terra::SpatRaster` with one allocation layer per time step.}
#'   \item{time_index}{Character vector of time labels.}
#'   \item{change_map_first_last}{A `terra::SpatRaster` indicating whether allocation
#'   changed between the first and last time step (`1` = changed, `0` = unchanged).}
#'   \item{persistence}{A `terra::SpatRaster` indicating whether each cell retained
#'   the same allocation across all time steps (`1` = persistent, `0` = changed at least once).}
#'   \item{polygons}{Optional list of `sf` polygon outputs by time.}
#'   \item{summaries}{Optional list of summary tables by time.}
#' }
#' @export
weighted_voronoi_time <- function(
    points_list,
    weight_col,
    boundary_sf,
    time_index = NULL,
    distance = c("euclidean", "geodesic"),
    geodesic_engine = c("multisource", "classic"),
    res = 20,
    resistance_list = NULL,
    dem_list = NULL,
    keep_polygons = FALSE,
    keep_summaries = TRUE,
    prepared = NULL,
    verbose = TRUE,
    ...
) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  
  distance <- match.arg(distance)
  geodesic_engine <- match.arg(geodesic_engine)
  
  if (!inherits(boundary_sf, "sf")) stop("boundary_sf must be an sf object.")
  if (!is.list(points_list) || length(points_list) < 1) {
    stop("points_list must be a non-empty list of sf POINT objects.")
  }
  if (!all(vapply(points_list, inherits, logical(1), "sf"))) {
    stop("All elements of points_list must be sf objects.")
  }
  
  n_t <- length(points_list)
  
  if (is.null(time_index)) {
    time_index <- names(points_list)
    if (is.null(time_index) || any(time_index == "")) {
      time_index <- paste0("t", seq_len(n_t))
    }
  }
  if (length(time_index) != n_t) {
    stop("time_index must have the same length as points_list.")
  }
  
  if (!is.null(resistance_list)) {
    if (!is.list(resistance_list)) stop("resistance_list must be a list if supplied.")
    if (!(length(resistance_list) %in% c(1, n_t))) {
      stop("resistance_list must have length 1 or the same length as points_list.")
    }
  }
  
  if (!is.null(dem_list)) {
    if (!is.list(dem_list)) stop("dem_list must be a list if supplied.")
    if (!(length(dem_list) %in% c(1, n_t))) {
      stop("dem_list must have length 1 or the same length as points_list.")
    }
  }
  
  dots <- list(...)
  dot1 <- function(name, default) {
    if (name %in% names(dots)) dots[[name]] else default
  }
  
  get_resistance_i <- function(i) {
    if (is.null(resistance_list)) return(NULL)
    if (length(resistance_list) == 1) resistance_list[[1]] else resistance_list[[i]]
  }
  
  get_dem_i <- function(i) {
    if (is.null(dem_list)) return(NULL)
    if (length(dem_list) == 1) dem_list[[1]] else dem_list[[i]]
  }
  
  allocation_list <- vector("list", n_t)
  polygon_list <- if (keep_polygons) vector("list", n_t) else NULL
  summary_list <- if (keep_summaries) vector("list", n_t) else NULL
  
  prepared_list <- NULL
  prepared_shared <- NULL
  
  if (distance == "geodesic") {
    static_surface <- (is.null(resistance_list) || length(resistance_list) == 1) &&
      (is.null(dem_list) || length(dem_list) == 1)
    
    if (static_surface) {
      prepared_shared <- prepare_geodesic_context(
        boundary_sf = boundary_sf,
        res = res,
        close_mask = dot1("close_mask", TRUE),
        close_iters = dot1("close_iters", 1),
        resistance_rast = get_resistance_i(1),
        dem_rast = get_dem_i(1),
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
    } else {
      prepared_list <- vector("list", n_t)
      
      for (i in seq_len(n_t)) {
        resistance_i <- get_resistance_i(i)
        dem_i <- get_dem_i(i)
        
        prepared_list[[i]] <- prepare_geodesic_context(
          boundary_sf = boundary_sf,
          res = res,
          close_mask = dot1("close_mask", TRUE),
          close_iters = dot1("close_iters", 1),
          resistance_rast = resistance_i,
          dem_rast = dem_i,
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
  }
  
  for (i in seq_len(n_t)) {
    if (verbose) {
      message(sprintf("Temporal tessellation %d / %d (%s)", i, n_t, time_index[i]))
    }
    
    points_i <- points_list[[i]]
    
    if (!weight_col %in% names(points_i)) {
      stop(sprintf("weight_col '%s' not found in points_list[[%d]].", weight_col, i))
    }
    
    resistance_i <- get_resistance_i(i)
    dem_i <- get_dem_i(i)
    prepared_i <- NULL
    if (distance == "geodesic") {
      prepared_i <- if (!is.null(prepared_shared)) prepared_shared else prepared_list[[i]]
    }
    
    out_i <- weighted_voronoi_domain(
      points_sf = points_i,
      weight_col = weight_col,
      boundary_sf = boundary_sf,
      distance = distance,
      geodesic_engine = geodesic_engine,
      res = res,
      resistance_rast = resistance_i,
      dem_rast = dem_i,
      prepared = prepared_i,
      verbose = FALSE,
      ...
    )
    
    allocation_i <- out_i$allocation
    names(allocation_i) <- time_index[i]
    allocation_list[[i]] <- allocation_i
    
    if (keep_polygons) {
      polygon_list[[i]] <- out_i$polygons
    }
    
    if (keep_summaries) {
      summary_list[[i]] <- out_i$summary
    }
  }
  
  allocation_stack <- terra::rast(allocation_list)
  names(allocation_stack) <- time_index
  
  # Change map: first vs last timestep
  if (terra::nlyr(allocation_stack) >= 2) {
    first_alloc <- allocation_stack[[1]]
    last_alloc  <- allocation_stack[[terra::nlyr(allocation_stack)]]
    
    change_map_first_last <- terra::ifel(first_alloc != last_alloc, 1, 0)
    change_map_first_last <- terra::ifel(
      is.na(first_alloc) | is.na(last_alloc),
      NA,
      change_map_first_last
    )
  } else {
    change_map_first_last <- allocation_stack[[1]]
    terra::values(change_map_first_last) <- 0
  }
  names(change_map_first_last) <- "change_map_first_last"
  
  # Persistence: same allocation across all time steps
  persistence <- terra::app(allocation_stack, fun = function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_real_)
    if (length(unique(x)) == 1) return(1)
    0
  })
  names(persistence) <- "persistence"
  
  out <- list(
    allocations = allocation_stack,
    time_index = time_index,
    change_map_first_last = change_map_first_last,
    persistence = persistence
  )
  
  if (keep_polygons) {
    names(polygon_list) <- time_index
    out$polygons <- polygon_list
  }
  
  if (keep_summaries) {
    names(summary_list) <- time_index
    out$summaries <- summary_list
  }
  
  out
}