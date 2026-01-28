build_summary <- function(polygons_sf,
                          allocation_rast,
                          points_sf_used,
                          weight_col,
                          weight_transform,
                          method) {
  
  stopifnot(inherits(polygons_sf, "sf"))
  stopifnot(inherits(allocation_rast, "SpatRaster"))
  stopifnot(inherits(points_sf_used, "sf"))
  
  # --- polygon areas by generator_id ---
  if (!"generator_id" %in% names(polygons_sf)) stop("polygons_sf must contain 'generator_id'")
  
  area_m2 <- as.numeric(sf::st_area(polygons_sf))
  gid <- polygons_sf$generator_id
  
  area_by_id <- tapply(area_m2, gid, sum, simplify = TRUE)
  total_area <- sum(area_m2)
  
  # --- polygon centroid per generator_id (use representative point to avoid centroids outside concave polys) ---
  rep_pts <- sf::st_point_on_surface(polygons_sf)
  coords <- sf::st_coordinates(rep_pts)
  # average reps if split features exist
  cx <- tapply(coords[,1], gid, mean, simplify = TRUE)
  cy <- tapply(coords[,2], gid, mean, simplify = TRUE)
  
  # --- raster cell counts per generator_id ---
  fr <- freq_noNA(allocation_rast)
  n_cells_by_id <- setNames(fr$count, fr$value)
  
  # --- weights ---
  w_raw <- as.numeric(points_sf_used[[weight_col]])
  w_used <- weight_transform(w_raw)
  
  # build table aligned by generator_id 1..n
  n <- nrow(points_sf_used)
  generator_id <- seq_len(n)
  
  out <- data.frame(
    generator_id = generator_id,
    point_row = generator_id,                 # row index within the USED points_sf
    weight_raw = w_raw,
    weight_used = w_used,
    area_m2 = as.numeric(area_by_id[as.character(generator_id)]),
    area_share = as.numeric(area_by_id[as.character(generator_id)]) / total_area,
    n_cells = as.integer(n_cells_by_id[as.character(generator_id)]),
    centroid_x = as.numeric(cx[as.character(generator_id)]),
    centroid_y = as.numeric(cy[as.character(generator_id)]),
    method = method,
    stringsAsFactors = FALSE
  )
  
  # Replace missing with 0 where appropriate (e.g. if something tiny got dropped before polygonise)
  out$area_m2[is.na(out$area_m2)] <- 0
  out$area_share[is.na(out$area_share)] <- 0
  out$n_cells[is.na(out$n_cells)] <- 0
  
  out
}


#----------------------------------------------------------------------------#


build_diagnostics <- function(boundary_sf,
                              polygons_sf,
                              points_sf_input,
                              points_sf_used,
                              allocation_rast,
                              method,
                              res,
                              unreachable_rast = NULL,
                              extra = list()) {
  
  stopifnot(inherits(boundary_sf, "sf"))
  stopifnot(inherits(polygons_sf, "sf"))
  stopifnot(inherits(allocation_rast, "SpatRaster"))
  
  boundary_area <- as.numeric(sf::st_area(sf::st_union(sf::st_make_valid(boundary_sf))))
  allocated_area <- sum(as.numeric(sf::st_area(polygons_sf)))
  coverage_rel_error <- if (boundary_area > 0) abs(boundary_area - allocated_area) / boundary_area else NA_real_
  
  # total allocated cells (inside mask)
  total_cells <- terra::global(!is.na(allocation_rast), "sum", na.rm = TRUE)[1,1]
  
  # unreachable info (geodesic only)
  n_unreachable <- NA_real_
  unreachable_fraction <- NA_real_
  if (!is.null(unreachable_rast)) {
    n_unreachable <- terra::global(unreachable_rast, "sum", na.rm = TRUE)[1,1]
    unreachable_fraction <- if (is.finite(total_cells) && total_cells > 0) n_unreachable / total_cells else NA_real_
  }
  
  diag <- list(
    method = method,
    res = res,
    n_points_input = if (!missing(points_sf_input) && !is.null(points_sf_input)) nrow(points_sf_input) else NA_integer_,
    n_points_used = if (!missing(points_sf_used) && !is.null(points_sf_used)) nrow(points_sf_used) else NA_integer_,
    n_points_dropped = if (!missing(points_sf_input) && !missing(points_sf_used) &&
                           !is.null(points_sf_input) && !is.null(points_sf_used))
      (nrow(points_sf_input) - nrow(points_sf_used)) else NA_integer_,
    boundary_area_m2 = boundary_area,
    allocated_area_m2 = allocated_area,
    coverage_rel_error = coverage_rel_error,
    total_cells = total_cells,
    n_unreachable_cells = n_unreachable,
    unreachable_fraction = unreachable_fraction
  )
  
  # attach any extra metadata (e.g. close_mask, close_iters, island settings, timing)
  if (length(extra)) {
    for (nm in names(extra)) diag[[nm]] <- extra[[nm]]
  }
  
  diag
}


#---------------------------------------------------------------------------------#


# Majority / modal value (ignores NA / non-finite)
mode1 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#---------------------------------------------------------------------------------#

# Force each generator to own at least the raster cell containing it
force_points_to_own_cell <- function(ID, points_sf) {
  stopifnot(inherits(ID, "SpatRaster"))
  pts_v <- terra::vect(points_sf)
  xy <- terra::crds(pts_v)
  cell <- terra::cellFromXY(ID, xy)
  
  v <- terra::values(ID)
  v[cell] <- seq_len(nrow(points_sf))
  terra::setValues(ID, v)
}

#---------------------------------------------------------------------------------#

freq_noNA <- function(x) {
  # Version-safe terra::freq() that drops NA rows if present
  fr <- terra::freq(x)
  if (is.null(fr) || nrow(fr) == 0) return(fr)
  
  # terra freq output column names vary; usually "value" and "count"
  if ("value" %in% names(fr)) {
    fr <- fr[is.finite(fr$value), , drop = FALSE]
  } else {
    # fallback: assume first column is the value
    fr <- fr[is.finite(fr[[1]]), , drop = FALSE]
    names(fr)[1] <- "value"
  }
  
  if (!"count" %in% names(fr)) {
    # fallback: assume second column is count
    if (ncol(fr) >= 2) names(fr)[2] <- "count"
  }
  
  fr
}

#---------------------------------------------------------------------------------#

# Remove disconnected "islands" so each generator has one connected region.
# Keeps the connected component that contains the generator point's cell.
# Reassigns removed cells by neighbor majority fill (and final nearest-neighbor fill if needed).
remove_islands <- function(ID, points_sf, min_cells = 5, max_iter_fill = 50) {
  stopifnot(inherits(ID, "SpatRaster"))
  
  ID2 <- ID
  
  pts_v <- terra::vect(points_sf)
  xy <- terra::crds(pts_v)
  pt_cells <- terra::cellFromXY(ID2, xy)
  
  labs <- sort(unique(terra::values(ID2)))
  labs <- labs[is.finite(labs)]
  
  for (lab in labs) {
    m <- (ID2 == lab)
    p <- terra::patches(m, directions = 8)
    
    fr <- freq_noNA(p)
    if (is.null(fr) || nrow(fr) == 0) next
    
    # Identify which patch contains the generator point for this label
    lab_index <- as.integer(lab)
    
    keep_patch <- NA_real_
    if (lab_index >= 1 && lab_index <= length(pt_cells) && is.finite(pt_cells[lab_index])) {
      keep_patch <- terra::values(p)[pt_cells[lab_index]]
    }
    
    # Fallback: keep largest patch
    if (!is.finite(keep_patch)) {
      fr <- fr[order(fr$count, decreasing = TRUE), ]
      keep_patch <- fr$value[1]
    }
    
    # Remove ALL other components (enforces "no islands")
    drop_patches <- fr$value[fr$value != keep_patch]
    
    # If you prefer to only remove small islands, use:
    # drop_patches <- fr$value[fr$value != keep_patch & fr$count < min_cells]
    
    if (length(drop_patches)) {
      ID2 <- terra::ifel(p %in% drop_patches, NA, ID2)
    }
  }
  
  # Fill NA holes by neighbor majority
  for (k in seq_len(max_iter_fill)) {
    n_na <- terra::global(is.na(ID2), "sum", na.rm = TRUE)[1, 1]
    if (is.na(n_na) || n_na == 0) break
    
    neigh <- terra::focal(
      ID2, w = matrix(1, 3, 3), fun = mode1,
      na.policy = "omit", fillvalue = NA
    )
    ID2 <- terra::ifel(is.na(ID2), neigh, ID2)
  }
  
  # Final safety net (version-safe): keep filling by local majority
  for (k in seq_len(200)) {
    n_na_final <- terra::global(is.na(ID2), "sum", na.rm = TRUE)[1, 1]
    if (is.na(n_na_final) || n_na_final == 0) break
    
    neigh <- terra::focal(
      ID2, w = matrix(1, 3, 3), fun = mode1,
      na.policy = "omit", fillvalue = NA
    )
    ID2 <- terra::ifel(is.na(ID2), neigh, ID2)
  }
  
  ID2
}

#---------------------------------------------------------------------------------#

#' Weighted Euclidean tessellation (core)
#'
#' Internal/core function used by [weighted_voronoi_domain()] to compute a weighted
#' Euclidean tessellation on a rasterised domain.
#'
#' @inheritParams weighted_voronoi_domain
#' @return A list containing polygon output (if requested), allocation raster, and weights.
#' @export

weighted_voronoi <- function(points_sf,
                             weight_col,
                             boundary = NULL,        # sf POLYGON/MULTIPOLYGON (optional)
                             template_rast = NULL,   # terra SpatRaster (optional)
                             res = NULL,             # numeric; only used if boundary is provided
                             weight_transform = function(w) w,
                             method = c("argmin", "partition"),
                             max_dist = NULL,
                             verbose = TRUE,
                             island_min_cells = 5,
                             island_fill_iter = 50) {
  
  method <- match.arg(method)
  
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  
  stopifnot(inherits(points_sf, "sf"))
  if (!weight_col %in% names(points_sf)) stop("weight_col not found in points_sf")
  
  boundary_use <- NULL  # keep a copy for later clipping
  
  # ---- choose raster template ----
  if (is.null(template_rast) && is.null(boundary)) {
    stop("Provide either boundary (sf polygon) or template_rast (terra SpatRaster).")
  }
  if (!is.null(template_rast) && !is.null(boundary)) {
    stop("Provide only one of boundary or template_rast, not both.")
  }
  
  if (!is.null(template_rast)) {
    r <- template_rast
  } else {
    if (is.null(res)) stop("If boundary is provided, res must be provided.")
    if (!inherits(boundary, "sf")) stop("boundary must be an sf polygon object.")
    
    if (sf::st_crs(points_sf) != sf::st_crs(boundary)) {
      boundary <- sf::st_transform(boundary, sf::st_crs(points_sf))
    }
    if (sf::st_is_longlat(points_sf)) {
      stop("Please use a projected CRS (meters), not lon/lat.")
    }
    
    boundary_use <- sf::st_make_valid(boundary)
    
    bnd_v <- terra::vect(boundary_use)
    r <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
    r <- terra::setValues(r, 1)   # important so mask() works on all terra builds
    r <- terra::crop(r, bnd_v)
    r <- terra::mask(r, bnd_v)
  }
  
  # ---- harmonise CRS if using template raster ----
  if (!is.null(template_rast)) {
    r_crs <- terra::crs(r, proj = TRUE)
    pts_crs <- sf::st_crs(points_sf)$wkt
    if (!identical(pts_crs, r_crs)) {
      points_sf <- sf::st_transform(points_sf, r_crs)
    }
    if (sf::st_is_longlat(points_sf)) stop("Template raster appears lon/lat; use a projected CRS.")
  }
  
  # ---- weights ----
  w_raw <- as.numeric(points_sf[[weight_col]])
  w <- weight_transform(w_raw)
  if (any(!is.finite(w)) || any(w <= 0)) stop("Transformed weights must be finite and > 0.")
  
  pts_v <- terra::vect(points_sf)
  
  # ---- build weighted distance minimum surface + winning id ----
  E  <- r; terra::values(E)  <- Inf
  ID <- r; terra::values(ID) <- NA_integer_
  
  for (i in seq_len(nrow(points_sf))) {
    if (verbose && (i %% 25 == 0 || i == 1 || i == nrow(points_sf))) {
      message(sprintf("Point %d / %d", i, nrow(points_sf)))
    }
    
    di <- terra::distance(r, pts_v[i])
    if (!is.null(max_dist)) di <- terra::clamp(di, upper = max_dist, values = TRUE)
    
    Ei <- di / w[i]
    upd <- Ei < E
    E  <- terra::ifel(upd, Ei, E)
    ID <- terra::ifel(upd, i,  ID)
  }
  
  if (method == "argmin") {
    return(list(allocation = ID, cost_surface = E, weights = w, weights_raw = w_raw))
  }
  
  # ---- Partition-based output + island removal ----
  ID <- force_points_to_own_cell(ID, points_sf)
  
  ID_clean <- remove_islands(
    ID,
    points_sf = points_sf,
    min_cells = island_min_cells,
    max_iter_fill = island_fill_iter
  )
  
  # Polygonise
  poly <- terra::as.polygons(ID, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  poly_sf <- sf::st_as_sf(poly)
  if (!is.null(boundary_use)) sf::st_crs(poly_sf) <- sf::st_crs(boundary_use)
  names(poly_sf)[names(poly_sf) == names(ID)] <- "generator_id"
  
  # Ensure CRS is carried through (terra->sf can drop it on some setups)
  sf::st_crs(poly_sf) <- sf::st_crs(boundary_use)
  
  
  # ---- NEW: hard clip to boundary for exact edge matching ----
  if (!is.null(boundary_use)) {
    poly_sf <- sf::st_make_valid(poly_sf)
    
    # Intersection can split features; clip then dissolve back to one per generator_id
    poly_sf <- sf::st_intersection(poly_sf, boundary_use)
    
    poly_sf <- poly_sf |>
      dplyr::group_by(generator_id) |>
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ .x[1]), do_union = TRUE) |>
      dplyr::ungroup()
  }
  
  # Attach weights + original point attributes
  poly_sf$weight <- w[poly_sf$generator_id]
  pts_df <- sf::st_drop_geometry(points_sf)
  pts_df$generator_id <- seq_len(nrow(points_sf))
  poly_sf <- merge(poly_sf, pts_df, by = "generator_id", all.x = TRUE)
  
  return(list(
    polygons = poly_sf,
    allocation = ID_clean,
    cost_surface = E,
    weights = w,
    weights_raw = w_raw
  ))
}

#--------------------------------------------------------------------------------#

#' Weighted geodesic tessellation (core)
#'
#' Computes a weighted tessellation using domain-constrained (geodesic) distances.
#' Distances are calculated as shortest-path distances through a rasterised domain mask.
#'
#' @inheritParams weighted_voronoi_domain
#' @return A list containing polygon output, allocation raster, and weights.
#' @export

weighted_voronoi_geodesic <- function(points_sf, weight_col, boundary_sf,
                                      res = 20,
                                      weight_transform = function(w) w,
                                      close_mask = TRUE,
                                      close_iters = 1,
                                      verbose = TRUE) {
  
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("raster", quietly = TRUE)) stop("Install raster.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  
  # --- CRS / validity ---
  if (sf::st_is_longlat(points_sf)) stop("Project to a metric CRS first.")
  if (sf::st_crs(points_sf) != sf::st_crs(boundary_sf)) {
    boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(points_sf))
  }
  boundary_sf <- sf::st_make_valid(boundary_sf)
  
  # Keep points that intersect (touch/inside) the boundary
  inside <- sf::st_intersects(points_sf, boundary_sf, sparse = FALSE)[, 1]
  if (any(!inside)) {
    message(sprintf("Dropping %d point(s) outside boundary.", sum(!inside)))
    points_sf <- points_sf[inside, , drop = FALSE]
  }
  if (nrow(points_sf) == 0) stop("No points remain inside the boundary.")
  
  # --- raster template + mask ---
  bnd_v <- terra::vect(boundary_sf)
  r <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
  r <- terra::setValues(r, 1)
  r <- terra::mask(r, bnd_v)
  
  # Optional: close tiny gaps in the raster mask (connectivity fix)
  if (close_mask) {
    m <- !is.na(r)
    for (k in seq_len(close_iters)) {
      m <- terra::focal(m, w = matrix(1, 3, 3), fun = max, na.policy = "omit", fillvalue = 0)
    }
    r <- terra::ifel(m == 1, 1, NA)
  }
  
  # --- weights ---
  w_raw <- as.numeric(points_sf[[weight_col]])
  w <- weight_transform(w_raw)
  if (any(!is.finite(w)) || any(w <= 0)) stop("Weights must be finite and > 0.")
  
  # --- conductance surface for gdistance ---
  r_in <- raster::raster(r)
  vv <- raster::values(r_in)
  vv[is.na(vv)] <- NA
  vv[!is.na(vv)] <- 1
  raster::values(r_in) <- vv
  
  tr <- gdistance::transition(r_in, function(x) 1, directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")
  
  pts_sp <- as(points_sf, "Spatial")
  
  # --- build weighted geodesic cost stack ---
  cost_list <- vector("list", nrow(points_sf))
  for (i in seq_len(nrow(points_sf))) {
    if (verbose && (i == 1 || i == nrow(points_sf) || i %% 10 == 0)) {
      message(sprintf("Geodesic distances: %d / %d", i, nrow(points_sf)))
    }
    
    d_i <- gdistance::accCost(tr, pts_sp[i, ])  # RasterLayer
    d_i <- terra::rast(d_i)                     # SpatRaster
    
    # Unreachable from this point => Inf (cannot win)
    d_i[is.na(d_i)] <- Inf
    
    cost_list[[i]] <- d_i / w[i]
  }
  
  S <- terra::rast(cost_list)
  
  # Cells unreachable from ALL points (all Inf)
  min_cost <- terra::app(S, fun = function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(Inf)
    min(vals)
  })
  all_unreachable <- is.infinite(min_cost)
  
  # Initial geodesic allocation
  ID <- terra::which.min(S)
  
  # Mark unreachable cells as NA (purely geodesic behaviour)
  ID <- terra::ifel(all_unreachable, NA, ID)
  
  # --- Fill unreachable cells locally (still "geodesic-consistent") ---
  # mode1 must exist in your environment
  for (k in seq_len(100)) {
    n_na <- terra::global(is.na(ID), "sum", na.rm = TRUE)[1, 1]
    if (is.na(n_na) || n_na == 0) break
    
    neigh <- terra::focal(ID, w = matrix(1, 3, 3), fun = mode1,
                          na.policy = "omit", fillvalue = NA)
    ID <- terra::ifel(is.na(ID), neigh, ID)
  }
  # Final safety net: iterative neighbor fill (version-safe)
  for (k in seq_len(200)) {
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
  
  # Enforce point-cell ownership + remove islands (your existing helpers)
  ID <- force_points_to_own_cell(ID, points_sf)
  ID <- remove_islands(ID, points_sf)
  
  # Polygonise
  poly <- terra::as.polygons(ID, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  poly_sf <- sf::st_as_sf(poly)
  sf::st_crs(poly_sf) <- sf::st_crs(boundary_sf)
  names(poly_sf)[names(poly_sf) == names(ID)] <- "generator_id"
  
  # Attach attributes
  poly_sf$weight <- w[poly_sf$generator_id]
  pts_df <- sf::st_drop_geometry(points_sf)
  pts_df$generator_id <- seq_len(nrow(points_sf))
  poly_sf <- merge(poly_sf, pts_df, by = "generator_id", all.x = TRUE)
  
  list(
    polygons = poly_sf,
    allocation = ID,
    weights = w,
    weights_raw = w_raw,
    unreachable = all_unreachable,
    points_used = points_sf
  )
}

#--------------------------------------------------------------------------------#

#' Weighted tessellation in a constrained polygon domain
#'
#' Creates a complete, connected tessellation of a polygonal domain using either
#' weighted Euclidean distance or weighted geodesic (domain-constrained) distance.
#' Weights are supplied as an attribute of generator points and can be transformed
#' by a user-defined function prior to allocation.
#'
#' @param points_sf An `sf` POINT object containing generator locations and attributes.
#' @param weight_col Character. Name of the weight column in `points_sf`.
#' @param boundary_sf An `sf` POLYGON/MULTIPOLYGON defining the domain.
#' @param res Numeric. Raster resolution in CRS units (e.g. metres).
#' @param weight_transform Function. Transforms weights before allocation. Must return
#'   finite, strictly positive values.
#' @param distance Character. One of `"euclidean"` or `"geodesic"`.
#' @param max_dist Optional numeric. Maximum Euclidean distance to consider (euclidean only).
#' @param island_min_cells Integer. Minimum patch size used in island removal.
#' @param island_fill_iter Integer. Maximum iterations for filling reassigned cells.
#' @param clip_to_boundary Logical. If `TRUE`, polygon output is intersected with the
#'   input boundary for exact edge matching (euclidean only).
#' @param close_mask Logical. If `TRUE`, applies a morphological closing to the raster
#'   mask (geodesic only).
#' @param close_iters Integer. Number of closing iterations (geodesic only).
#' @param verbose Logical. If `TRUE`, prints progress.
#'
#' @return A list with elements including:
#' \describe{
#'   \item{polygons}{An `sf` object with one polygon per generator.}
#'   \item{allocation}{A `terra::SpatRaster` assigning each cell to a generator.}
#'   \item{summary}{A generator-level summary table.}
#'   \item{diagnostics}{A list of diagnostic metrics and settings.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' crs_use <- 32636
#' boundary_sf <- st_sf(
#'   geometry = st_sfc(st_polygon(list(rbind(
#'     c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
#'   )))),
#'   crs = crs_use
#' )
#' points_sf <- st_sf(
#'   population = c(50, 200, 1000),
#'   geometry = st_sfc(st_point(c(200,200)), st_point(c(800,250)), st_point(c(500,500))),
#'   crs = crs_use
#' )
#' out <- weighted_voronoi_domain(points_sf, "population", boundary_sf,
#'   res = 20, weight_transform = log10, distance = "euclidean", verbose = FALSE
#' )
#' }
#' @export

weighted_voronoi_domain <- function(points_sf,
                                    weight_col,
                                    boundary_sf,
                                    res = 20,
                                    weight_transform = function(w) w,
                                    distance = c("euclidean", "geodesic"),
                                    # Euclidean options
                                    max_dist = NULL,
                                    island_min_cells = 5,
                                    island_fill_iter = 50,
                                    clip_to_boundary = TRUE,
                                    # Geodesic options
                                    close_mask = TRUE,
                                    close_iters = 1,
                                    # general
                                    verbose = TRUE) {
  
  distance <- match.arg(distance)
  
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  
  stopifnot(inherits(points_sf, "sf"))
  stopifnot(inherits(boundary_sf, "sf"))
  if (!weight_col %in% names(points_sf)) stop("weight_col not found in points_sf")
  if (sf::st_is_longlat(points_sf)) stop("Please use a projected CRS (meters), not lon/lat.")
  
  points_input <- points_sf
  
  # Harmonise CRS
  boundary_sf <- sf::st_make_valid(boundary_sf)
  if (sf::st_crs(points_sf) != sf::st_crs(boundary_sf)) {
    boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(points_sf))
  }
  
  if (distance == "euclidean") {
    
    out <- weighted_voronoi(
      points_sf = points_sf,
      weight_col = weight_col,
      boundary = boundary_sf,
      res = res,
      weight_transform = weight_transform,
      method = "partition",
      max_dist = max_dist,
      verbose = verbose,
      island_min_cells = island_min_cells,
      island_fill_iter = island_fill_iter
    )
    
    # Optional exact clip (for perfect edge match in comparisons)
    if (clip_to_boundary && !is.null(boundary_sf)) {
      poly_sf <- sf::st_make_valid(out$polygons)
      b <- sf::st_make_valid(boundary_sf)
      
      poly_sf <- sf::st_intersection(poly_sf, b)
      
      # If intersection produced nothing (rare), stop early with a clear message
      if (nrow(poly_sf) == 0) stop("Clipping produced 0 polygons. Check boundary CRS/validity.")
      
      gids <- sort(unique(poly_sf$generator_id))
      
      geoms <- lapply(gids, function(g) {
        u <- sf::st_union(sf::st_geometry(poly_sf[poly_sf$generator_id == g, , drop = FALSE]))
        # st_union returns an sfc of length 1; extract the single sfg:
        u[[1]]
      })
      
      out$polygons <- sf::st_sf(
        generator_id = gids,
        geometry = sf::st_sfc(geoms, crs = sf::st_crs(poly_sf))
      )
    }
    
    
    summary_tbl <- build_summary(
      polygons_sf = out$polygons,
      allocation_rast = out$allocation,
      points_sf_used = points_sf,
      weight_col = weight_col,
      weight_transform = weight_transform,
      method = "euclidean"
    )
    
    diagnostics <- build_diagnostics(
      boundary_sf = boundary_sf,
      polygons_sf = out$polygons,
      points_sf_input = points_input,
      points_sf_used = points_sf,
      allocation_rast = out$allocation,
      method = "euclidean",
      res = res,
      unreachable_rast = NULL,
      extra = list(
        island_min_cells = island_min_cells,
        island_fill_iter = island_fill_iter,
        clip_to_boundary = clip_to_boundary
      )
    )
    
    return(list(
      method = "euclidean",
      polygons = out$polygons,
      allocation = out$allocation,
      cost_surface = out$cost_surface,
      weights = out$weights,
      weights_raw = out$weights_raw,
      summary = summary_tbl,
      diagnostics = diagnostics
    ))
  }
  
  # distance == "geodesic"
  out <- weighted_voronoi_geodesic(
    points_sf = points_sf,
    weight_col = weight_col,
    boundary_sf = boundary_sf,
    res = res,
    weight_transform = weight_transform,
    close_mask = close_mask,
    close_iters = close_iters,
    verbose = verbose
  )
  
  # IMPORTANT: geodesic function may drop outside points internally; treat that as "used"
  # If your geodesic function currently drops points, consider returning points_used too.
  points_used <- points_sf
  if (!is.null(out$points_used)) points_used <- out$points_used
  
  summary_tbl <- build_summary(
    polygons_sf = out$polygons,
    allocation_rast = out$allocation,
    points_sf_used = points_used,
    weight_col = weight_col,
    weight_transform = weight_transform,
    method = "geodesic"
  )
  
  diagnostics <- build_diagnostics(
    boundary_sf = boundary_sf,
    polygons_sf = out$polygons,
    points_sf_input = points_input,
    points_sf_used = points_used,
    allocation_rast = out$allocation,
    method = "geodesic",
    res = res,
    unreachable_rast = out$unreachable,
    extra = list(
      close_mask = close_mask,
      close_iters = close_iters
    )
  )
  
  list(
    method = "geodesic",
    polygons = out$polygons,
    allocation = out$allocation,
    weights = out$weights,
    weights_raw = out$weights_raw,
    unreachable = out$unreachable,
    summary = summary_tbl,
    diagnostics = diagnostics
  )
}
