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
  
  e <- terra::ext(bnd_v)
  e[1] <- e[1] - res
  e[2] <- e[2] + res
  e[3] <- e[3] - res
  e[4] <- e[4] + res
  
  r <- terra::rast(
    ext = e,
    resolution = res,
    crs = terra::crs(bnd_v)
  )
  r <- terra::setValues(r, 1)
  r <- terra::mask(r, bnd_v)
  
  if (close_mask) {
    m <- !is.na(r)
    
    # dilation
    for (k in seq_len(close_iters)) {
      m <- terra::focal(
        m,
        w = matrix(1, 3, 3),
        fun = max,
        na.policy = "omit",
        fillvalue = 0
      )
    }
    
    # erosion
    for (k in seq_len(close_iters)) {
      m <- terra::focal(
        m,
        w = matrix(1, 3, 3),
        fun = min,
        na.policy = "omit",
        fillvalue = 0
      )
    }
    
    r <- terra::ifel(m == 1, 1, NA)
    r <- terra::mask(r, bnd_v)
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
                                         template_rast = NULL,
                                         verbose = TRUE) {
  stopifnot(inherits(points_sf, "sf"))
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  
  weight_model <- match.arg(weight_model, c("multiplicative", "power", "additive"))
  
  if (length(weight_vec) != nrow(points_sf)) {
    stop("weight_vec must have length equal to nrow(points_sf).")
  }
  if (is.null(template_rast) || !inherits(template_rast, "SpatRaster")) {
    stop("template_rast must be supplied as a terra SpatRaster.")
  }
  
  pts_sp <- methods::as(points_sf, "Spatial")
  cost_list <- vector("list", nrow(points_sf))
  
  inside_domain <- !is.na(terra::values(template_rast, mat = FALSE))
  
  for (i in seq_len(nrow(points_sf))) {
    if (verbose && (i == 1 || i == nrow(points_sf) || i %% 10 == 0)) {
      message(sprintf("Geodesic distances: %d / %d", i, nrow(points_sf)))
    }
    
    d_i <- gdistance::accCost(tr, pts_sp[i, ])
    d_i <- terra::rast(d_i)
    
    d_vals <- terra::values(d_i, mat = FALSE)
    
    # Inside the domain but unreachable from this source -> Inf
    d_vals[inside_domain & is.na(d_vals)] <- Inf
    
    # Outside the true domain mask must stay NA
    d_vals[!inside_domain] <- NA_real_
    
    d_i <- terra::setValues(d_i, d_vals)
    
    cost_i <- .effective_cost(
      d = d_i,
      w = weight_vec[i],
      weight_model = weight_model,
      weight_power = weight_power
    )
    
    # Safety: keep outside-domain cells as NA after weighting too
    cost_vals <- terra::values(cost_i, mat = FALSE)
    cost_vals[!inside_domain] <- NA_real_
    cost_i <- terra::setValues(cost_i, cost_vals)
    
    cost_list[[i]] <- cost_i
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

.fill_na_by_neighbors <- function(ID, domain_mask, max_iter = 300) {
  stopifnot(inherits(ID, "SpatRaster"))
  stopifnot(inherits(domain_mask, "SpatRaster"))
  
  mask_vals <- terra::values(domain_mask, mat = FALSE)
  
  for (k in seq_len(max_iter)) {
    vals <- terra::values(ID, mat = FALSE)
    n_na <- sum(is.na(vals) & mask_vals == 1, na.rm = TRUE)
    if (n_na == 0) break
    
    neigh <- terra::focal(
      ID,
      w = matrix(1, 3, 3),
      fun = mode1,
      na.policy = "omit",
      fillvalue = NA
    )
    
    neigh_vals <- terra::values(neigh, mat = FALSE)
    fill_idx <- which(is.na(vals) & mask_vals == 1 & !is.na(neigh_vals))
    
    if (!length(fill_idx)) break
    
    vals[fill_idx] <- neigh_vals[fill_idx]
    vals[mask_vals != 1] <- NA
    ID <- terra::setValues(ID, vals)
  }
  
  ID
}

.postprocess_allocation <- function(ID,
                                    points_sf,
                                    domain_mask,
                                    island_min_cells = 5,
                                    island_fill_iter = 50,
                                    fill_iter = 300) {
  stopifnot(inherits(ID, "SpatRaster"))
  stopifnot(inherits(points_sf, "sf"))
  stopifnot(inherits(domain_mask, "SpatRaster"))
  
  mask_vals <- terra::values(domain_mask, mat = FALSE)
  
  # Enforce mask before anything else
  vals <- terra::values(ID, mat = FALSE)
  vals[mask_vals != 1] <- NA
  ID <- terra::setValues(ID, vals)
  
  # Fill internal gaps only inside the domain
  ID <- .fill_na_by_neighbors(ID, domain_mask = domain_mask, max_iter = fill_iter)
  
  # Ensure each point owns its containing cell
  ID <- force_points_to_own_cell(ID, points_sf)
  
  # Re-enforce mask
  vals <- terra::values(ID, mat = FALSE)
  vals[mask_vals != 1] <- NA
  ID <- terra::setValues(ID, vals)
  
  # Remove disconnected islands
  ID <- remove_islands(
    ID,
    points_sf = points_sf,
    domain_mask = domain_mask,
    min_cells = island_min_cells,
    max_iter_fill = island_fill_iter
  )
  
  # Re-enforce mask one final time
  vals <- terra::values(ID, mat = FALSE)
  vals[mask_vals != 1] <- NA
  ID <- terra::setValues(ID, vals)
  
  ID
}