# Internal helpers for anisotropic geodesic tessellation
# No exports

.align_dem_to_mask <- function(dem_rast, mask_r) {
  stopifnot(inherits(dem_rast, "SpatRaster"), inherits(mask_r, "SpatRaster"))
  
  if (!terra::same.crs(dem_rast, mask_r)) {
    dem_rast <- terra::project(dem_rast, mask_r)
  }
  if (!terra::compareGeom(dem_rast, mask_r, stopOnError = FALSE)) {
    dem_rast <- terra::resample(dem_rast, mask_r, method = "bilinear")
  }
  
  terra::mask(dem_rast, mask_r)
}

.cell_xy_matrix <- function(r) {
  stopifnot(inherits(r, "SpatRaster"))
  terra::xyFromCell(r, seq_len(terra::ncell(r)))
}

.build_neighbor_table <- function(r, directions = 8) {
  stopifnot(inherits(r, "SpatRaster"))
  
  # terra::adjacent returns from/to cell indices
  adj <- terra::adjacent(
    x = r,
    cells = seq_len(terra::ncell(r)),
    directions = directions,
    pairs = TRUE,
    include = FALSE
  )
  
  as.data.frame(adj)
}

.compute_directional_conductance <- function(
    from_xy, to_xy,
    z_from, z_to,
    v0_kmh = 6,
    a = 3.5,
    b = 0.05,
    min_speed_kmh = 0.25,
    uphill_factor = 1,
    downhill_factor = 1) {
  
  dx <- to_xy[, 1] - from_xy[, 1]
  dy <- to_xy[, 2] - from_xy[, 2]
  dxy <- sqrt(dx^2 + dy^2)
  
  # guard against zero-length edges
  bad_d <- !is.finite(dxy) | dxy <= 0
  dxy[bad_d] <- NA_real_
  
  dz <- z_to - z_from
  slope <- dz / dxy   # signed slope
  
  # Signed Tobler-style speed:
  # fastest on mild downhill, slower on steep uphill/downhill
  v_kmh <- v0_kmh * exp(-a * abs(slope + b))
  
  # Optional extra directional asymmetry controls
  uphill_idx <- is.finite(slope) & slope > 0
  downhill_idx <- is.finite(slope) & slope < 0
  
  if (any(uphill_idx)) {
    v_kmh[uphill_idx] <- v_kmh[uphill_idx] / uphill_factor
  }
  if (any(downhill_idx)) {
    v_kmh[downhill_idx] <- v_kmh[downhill_idx] * downhill_factor
  }
  
  # Clamp minimum speed
  v_kmh[!is.finite(v_kmh)] <- NA_real_
  v_kmh <- pmax(v_kmh, min_speed_kmh, na.rm = FALSE)
  
  # Convert to m/s
  v_ms <- v_kmh * (1000 / 3600)
  
  # Conductance as ease of travel across the edge
  conductance <- v_ms / dxy
  
  # clean invalid values
  conductance[!is.finite(conductance)] <- 0
  conductance[conductance < 0] <- 0
  
  conductance
}

.build_terrain_anisotropic_transition <- function(dem_rast,
                                                  mask_r,
                                                  uphill_factor = 1,
                                                  downhill_factor = 1,
                                                  v0_kmh = 6,
                                                  a = 3.5,
                                                  b = 0.05,
                                                  min_speed_kmh = 0.25) {
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("raster", quietly = TRUE)) stop("Install raster.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix.")
  
  stopifnot(inherits(dem_rast, "SpatRaster"), inherits(mask_r, "SpatRaster"))
  
  if (!is.finite(uphill_factor) || uphill_factor <= 0) {
    stop("uphill_factor must be finite and > 0.")
  }
  if (!is.finite(downhill_factor) || downhill_factor <= 0) {
    stop("downhill_factor must be finite and > 0.")
  }
  
  dem <- .align_dem_to_mask(dem_rast, mask_r)
  
  # Valid cells only within domain
  mask_vals <- terra::values(mask_r, mat = FALSE)
  valid_cells <- which(!is.na(mask_vals))
  
  # Build all neighbor pairs over full raster
  adj <- .build_neighbor_table(mask_r, directions = 8)
  names(adj) <- c("from", "to")
  
  # Keep only pairs where both cells are inside domain
  keep <- adj$from %in% valid_cells & adj$to %in% valid_cells
  adj <- adj[keep, , drop = FALSE]
  
  # Coordinates and elevations
  xy <- .cell_xy_matrix(mask_r)
  z <- terra::values(dem, mat = FALSE)
  
  from_xy <- xy[adj$from, , drop = FALSE]
  to_xy   <- xy[adj$to, , drop = FALSE]
  z_from  <- z[adj$from]
  z_to    <- z[adj$to]
  
  # Drop any pairs with missing elevation
  ok <- is.finite(z_from) & is.finite(z_to)
  adj <- adj[ok, , drop = FALSE]
  from_xy <- from_xy[ok, , drop = FALSE]
  to_xy   <- to_xy[ok, , drop = FALSE]
  z_from  <- z_from[ok]
  z_to    <- z_to[ok]
  
  conductance <- .compute_directional_conductance(
    from_xy = from_xy,
    to_xy = to_xy,
    z_from = z_from,
    z_to = z_to,
    v0_kmh = v0_kmh,
    a = a,
    b = b,
    min_speed_kmh = min_speed_kmh,
    uphill_factor = uphill_factor,
    downhill_factor = downhill_factor
  )
  
  # Safety cleanup
  conductance[!is.finite(conductance)] <- 0
  conductance[conductance < 0] <- 0
  
  n <- terra::ncell(mask_r)
  
  M <- Matrix::sparseMatrix(
    i = adj$from,
    j = adj$to,
    x = conductance,
    dims = c(n, n)
  )
  
  r_mask <- raster::raster(mask_r)
  
  # Build a valid TransitionLayer template, then replace its matrix
  tr <- gdistance::transition(
    x = r_mask,
    transitionFunction = function(x) 1,
    directions = 8,
    symm = FALSE
  )
  
  gdistance::transitionMatrix(tr) <- M
  tr@matrixValues <- "conductance"
  
  tr <- gdistance::geoCorrection(tr, type = "c")
  
  tr
}