# Internal helpers for geodesic tessellation
# No exports

.build_domain_mask <- function(boundary_sf,
                               res,
                               close_mask = TRUE,
                               close_iters = 1) {
  stopifnot(inherits(boundary_sf, "sf"))
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  
  boundary_sf <- sf::st_make_valid(boundary_sf)
  bnd_v <- terra::vect(boundary_sf)
  
  r <- terra::rast(
    ext = terra::ext(bnd_v),
    resolution = res,
    crs = terra::crs(bnd_v)
  )
  r <- terra::setValues(r, 1)
  r <- terra::mask(r, bnd_v)
  
  if (close_mask) {
    m <- !is.na(r)
    for (k in seq_len(close_iters)) {
      m <- terra::focal(
        m,
        w = matrix(1, 3, 3),
        fun = max,
        na.policy = "omit",
        fillvalue = 0
      )
    }
    r <- terra::ifel(m == 1, 1, NA)
  }
  
  r
}

.build_isotropic_transition <- function(resistance_rast) {
  stopifnot(inherits(resistance_rast, "SpatRaster"))
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("raster", quietly = TRUE)) stop("Install raster.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  
  cond <- 1 / resistance_rast
  cond[!is.finite(cond)] <- 0
  
  r_in <- raster::raster(cond)
  
  tr <- gdistance::transition(
    r_in,
    transitionFunction = mean,
    directions = 8
  )
  tr <- gdistance::geoCorrection(tr, type = "c")
  
  tr
}

.compute_geodesic_cost_stack <- function(tr,
                                         points_sf,
                                         weight_vec,
                                         weight_model = "multiplicative",
                                         weight_power = 1,
                                         verbose = TRUE) {
  stopifnot(inherits(points_sf, "sf"))
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  
  weight_model <- match.arg(weight_model, c("multiplicative", "power", "additive"))
  
  if (length(weight_vec) != nrow(points_sf)) {
    stop("weight_vec must have length equal to nrow(points_sf).")
  }
  
  pts_sp <- methods::as(points_sf, "Spatial")
  cost_list <- vector("list", nrow(points_sf))
  
  for (i in seq_len(nrow(points_sf))) {
    if (verbose && (i == 1 || i == nrow(points_sf) || i %% 10 == 0)) {
      message(sprintf("Geodesic distances: %d / %d", i, nrow(points_sf)))
    }
    
    d_i <- gdistance::accCost(tr, pts_sp[i, ])
    d_i <- terra::rast(d_i)
    
    # Unreachable from this source cannot win
    d_i[is.na(d_i)] <- Inf
    
    cost_list[[i]] <- .effective_cost(
      d = d_i,
      w = weight_vec[i],
      weight_model = weight_model,
      weight_power = weight_power
    )
  }
  
  terra::rast(cost_list)
}

.detect_all_unreachable <- function(cost_stack) {
  stopifnot(inherits(cost_stack, "SpatRaster"))
  
  min_cost <- terra::app(cost_stack, fun = function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(Inf)
    min(vals)
  })
  
  is.infinite(min_cost)
}

.assign_from_cost_stack <- function(cost_stack, all_unreachable = NULL) {
  stopifnot(inherits(cost_stack, "SpatRaster"))
  
  ID <- terra::which.min(cost_stack)
  
  if (!is.null(all_unreachable)) {
    ID <- terra::ifel(all_unreachable, NA, ID)
  }
  
  ID
}

.fill_na_by_neighbors <- function(ID, max_iter = 300) {
  stopifnot(inherits(ID, "SpatRaster"))
  
  for (k in seq_len(max_iter)) {
    n_na <- terra::global(is.na(ID), "sum", na.rm = TRUE)[1, 1]
    if (is.na(n_na) || n_na == 0) break
    
    neigh <- terra::focal(
      ID,
      w = matrix(1, 3, 3),
      fun = mode1,
      na.policy = "omit",
      fillvalue = NA
    )
    
    ID <- terra::ifel(is.na(ID), neigh, ID)
  }
  
  ID
}

.postprocess_allocation <- function(ID,
                                    points_sf,
                                    island_min_cells = 5,
                                    island_fill_iter = 50,
                                    fill_iter = 300) {
  stopifnot(inherits(ID, "SpatRaster"))
  stopifnot(inherits(points_sf, "sf"))
  
  # Fill local gaps first
  ID <- .fill_na_by_neighbors(ID, max_iter = fill_iter)
  
  # Ensure each point owns its containing cell
  ID <- force_points_to_own_cell(ID, points_sf)
  
  # Remove disconnected islands / small patches
  ID <- remove_islands(
    ID,
    points_sf = points_sf,
    min_cells = island_min_cells,
    max_iter_fill = island_fill_iter
  )
  
  ID
}