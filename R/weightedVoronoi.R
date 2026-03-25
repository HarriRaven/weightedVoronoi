.effective_cost <- function(d, w, weight_model = c("multiplicative", "power", "additive"), weight_power = 1) {
  weight_model <- match.arg(weight_model)
  
  if (!is.finite(w) || w <= 0) stop("Weights must be finite and > 0.")
  if (weight_model == "power" && (!is.finite(weight_power) || weight_power <= 0)) {
    stop("weight_power must be finite and > 0 when weight_model='power'.")
  }
  
  if (weight_model == "multiplicative") {
    return(d / w)
  }
  if (weight_model == "power") {
    return((d ^ weight_power) / w)
  }
  # additive
  return(d + (1 / w))
}

.prep_resistance <- function(mask_r,
                             resistance_rast = NULL,
                             dem_rast = NULL,
                             use_tobler = TRUE,
                             tobler_v0_kmh = 6,
                             tobler_a = 3.5,
                             tobler_b = 0.05,
                             min_speed_kmh = 0.25) {
  # Default: uniform resistance (=1) across the domain mask
  resistance <- mask_r
  terra::values(resistance) <- 1
  
  # 1) User-supplied resistance takes priority
  if (!is.null(resistance_rast)) {
    if (!inherits(resistance_rast, "SpatRaster")) stop("resistance_rast must be a terra SpatRaster.")
    rr <- resistance_rast
    
    # Align CRS and grid to the mask raster
    if (!identical(terra::crs(rr, proj = TRUE), terra::crs(mask_r, proj = TRUE))) {
      rr <- terra::project(rr, mask_r)
    }
    rr <- terra::resample(rr, mask_r, method = "bilinear")
    rr <- terra::mask(rr, mask_r)
    
    # Validate positivity
    bad <- terra::global(rr <= 0, "sum", na.rm = TRUE)[1, 1]
    if (is.finite(bad) && bad > 0) stop("resistance_rast must be strictly > 0 within the domain.")
    return(rr)
  }
  
  # 2) Otherwise DEM -> Tobler (if DEM is supplied)
  if (!is.null(dem_rast)) {
    if (!inherits(dem_rast, "SpatRaster")) stop("dem_rast must be a terra SpatRaster.")
    if (!use_tobler) stop("dem_rast provided but use_tobler=FALSE; supply resistance_rast instead.")
    
    return(.tobler_resistance_from_dem(
      dem_rast = dem_rast,
      mask_r = mask_r,
      v0_kmh = tobler_v0_kmh,
      a = tobler_a,
      b = tobler_b,
      min_speed_kmh = min_speed_kmh
    ))
  }
  
  resistance
}

#' Compose a resistance surface from multiple raster layers
#'
#' Aligns one or more `terra::SpatRaster` layers to a common template (CRS, resolution,
#' extent) and combines them into a single resistance raster using a specified rule.
#' Intended for building `resistance_rast` inputs for geodesic tessellations.
#'
#' @param ... One or more `terra::SpatRaster` layers. All layers must represent strictly
#'   positive resistance values.
#' @param template Optional `terra::SpatRaster` defining the target grid (CRS, resolution,
#'   extent). Defaults to the first layer.
#' @param mask Optional mask applied after alignment. Can be a `terra::SpatRaster`,
#'   `terra::SpatVector`, or `sf` polygon.
#' @param method How to combine layers: `"multiply"` (default), `"add"`, or `"max"`.
#' @param resample_method Resampling method used when aligning layers: `"bilinear"` for
#'   continuous resistance, `"near"` for categorical rasters.
#' @param na_policy How to treat NA values: `"propagate"` makes any NA propagate to the
#'   output; `"ignore"` drops NAs (neutral element for the chosen method).
#'
#' @return A `terra::SpatRaster` resistance surface on the template grid.
#' @export
#'
#' @examples
#' \dontrun{
#' R <- compose_resistance(slope_resistance, landcover_resistance, method = "multiply")
#' R <- add_barriers(R, rivers_sf, permeability = "semi", cost_multiplier = 20, width = 30)
#' out <- weighted_voronoi_domain(points_sf, "w", boundary_sf, distance = "geodesic",
#'                                resistance_rast = R)
#' }
compose_resistance <- function(...,
                               template = NULL,
                               mask = NULL,
                               method = c("multiply", "add", "max"),
                               resample_method = c("bilinear", "near"),
                               na_policy = c("propagate", "ignore")) {
  
  method <- match.arg(method)
  resample_method <- match.arg(resample_method)
  na_policy <- match.arg(na_policy)
  
  layers <- list(...)
  layers <- layers[!vapply(layers, is.null, logical(1))]
  if (length(layers) < 1) stop("Provide at least one SpatRaster layer.")
  if (!all(vapply(layers, inherits, logical(1), "SpatRaster"))) {
    stop("All inputs to compose_resistance() must be terra SpatRaster objects.")
  }
  
  if (is.null(template)) template <- layers[[1]]
  
  align1 <- function(x) {
    if (!identical(terra::crs(x, proj = TRUE), terra::crs(template, proj = TRUE))) {
      x <- terra::project(x, template)
    }
    x <- terra::resample(x, template, method = resample_method)
    
    if (!is.null(mask)) {
      # mask can be SpatRaster or SpatVector/sf polygon
      x <- terra::mask(x, mask)
    }
    x
  }
  
  aligned <- lapply(layers, align1)
  S <- terra::rast(aligned)
  
  fun <- switch(
    method,
    multiply = function(...) {
      v <- c(...)
      if (na_policy == "propagate" && anyNA(v)) return(NA_real_)
      if (na_policy == "ignore") v <- v[!is.na(v)]
      if (!length(v)) return(NA_real_)
      prod(v)
    },
    add = function(...) {
      v <- c(...)
      if (na_policy == "propagate" && anyNA(v)) return(NA_real_)
      if (na_policy == "ignore") v <- v[!is.na(v)]
      if (!length(v)) return(NA_real_)
      sum(v)
    },
    max = function(...) {
      v <- c(...)
      if (na_policy == "propagate" && anyNA(v)) return(NA_real_)
      if (na_policy == "ignore") v <- v[!is.na(v)]
      if (!length(v)) return(NA_real_)
      max(v)
    }
  )
  
  out <- terra::app(S, fun = fun)
  
  # Validate and enforce strictly positive
  bad <- terra::global(out <= 0, "sum", na.rm = TRUE)[1, 1]
  if (is.finite(bad) && bad > 0) {
    # either stop or clamp; I'd default to stop for scientific correctness
    stop("Composed resistance has values <= 0. Ensure all layers are strictly > 0.")
  }
  
  out
}

.barrier_mask_raster <- function(barriers,
                                 template_rast,
                                 width = 0) {
  if (!inherits(template_rast, "SpatRaster")) stop("template_rast must be a terra SpatRaster.")
  
  # If barriers is already a raster, align it and convert to mask
  if (inherits(barriers, "SpatRaster")) {
    b <- barriers
    if (!identical(terra::crs(b, proj = TRUE), terra::crs(template_rast, proj = TRUE))) {
      b <- terra::project(b, template_rast)
    }
    b <- terra::resample(b, template_rast, method = "near")
    # Anything >0 treated as barrier
    m <- terra::ifel(b > 0, 1, 0)
    return(m)
  }
  
  # Convert sf -> SpatVector
  if (inherits(barriers, "sf")) {
    barriers <- terra::vect(barriers)
  }
  if (!inherits(barriers, "SpatVector")) {
    stop("barriers must be a SpatRaster, sf object, or terra SpatVector.")
  }
  
  # Reproject vector to template CRS if needed
  bcrs <- terra::crs(barriers, proj = TRUE)
  tcrs <- terra::crs(template_rast, proj = TRUE)
  if (!identical(bcrs, tcrs)) {
    barriers <- terra::project(barriers, tcrs)
  }
  
  # Optional buffering of line features to ensure they hit raster cells
  if (!is.null(width) && is.finite(width) && width > 0) {
    barriers <- terra::buffer(barriers, width = width)
  }
  
  # Rasterise to mask (1 where barrier exists)
  m <- terra::rasterize(barriers, template_rast, field = 1, background = 0, touches = TRUE)
  m <- terra::ifel(m >= 1, 1, 0)
  m
}

#' Add barriers to a resistance surface
#'
#' Modifies a resistance raster by applying semi-permeable or impermeable barriers
#' provided as a raster mask (values > 0 treated as barrier) or as vector features
#' (sf / SpatVector) rasterised onto the resistance grid.
#'
#' @param resistance terra::SpatRaster of strictly positive movement resistance.
#' @param barriers A terra::SpatRaster mask (values > 0 treated as barrier), or
#'   an sf/SpatVector LINESTRING/POLYGON object.
#' @param permeability One of "semi", "impermeable", "permeable".
#' @param cost_multiplier Numeric > 0. Multiplier applied where barrier present
#'   (for "semi" and "permeable").
#' @param width Buffer distance (CRS units) applied to vector barriers before rasterising.
#' @return A terra::SpatRaster resistance surface with barrier effects applied.
#' @export

add_barriers <- function(resistance,
                         barriers,
                         permeability = c("semi", "impermeable", "permeable"),
                         cost_multiplier = 10,
                         width = 0) {
  
  if (!inherits(resistance, "SpatRaster")) stop("resistance must be a terra SpatRaster.")
  permeability <- match.arg(permeability)
  
  if (!is.finite(cost_multiplier) || cost_multiplier <= 0) {
    stop("cost_multiplier must be finite and > 0.")
  }
  
  mask <- .barrier_mask_raster(barriers, template_rast = resistance, width = width)
  
  out <- resistance
  
  if (permeability == "permeable") {
    # small friction increase (still useful to let users mark features)
    out <- terra::ifel(mask == 1, out * max(1, cost_multiplier), out)
  } else if (permeability == "semi") {
    out <- terra::ifel(mask == 1, out * cost_multiplier, out)
  } else {
    # impermeable: encode as infinite resistance
    out <- terra::ifel(mask == 1, Inf, out)
  }
  
  out
}

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
  rep_geom <- sf::st_point_on_surface(sf::st_geometry(polygons_sf))
  coords <- sf::st_coordinates(rep_geom)
  # average reps if split features exist
  cx <- tapply(coords[,1], gid, mean, simplify = TRUE)
  cy <- tapply(coords[,2], gid, mean, simplify = TRUE)
  
  # --- raster cell counts per generator_id ---
  fr <- freq_noNA(allocation_rast)
  n_cells_by_id <- stats::setNames(fr$count, fr$value)
  
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
remove_islands <- function(ID, points_sf, domain_mask, min_cells = 5, max_iter_fill = 50) {
  stopifnot(inherits(ID, "SpatRaster"))
  stopifnot(inherits(domain_mask, "SpatRaster"))
  
  ID2 <- ID
  mask_vals <- terra::values(domain_mask, mat = FALSE)
  
  pts_v <- terra::vect(points_sf)
  xy <- terra::crds(pts_v)
  pt_cells <- terra::cellFromXY(ID2, xy)
  
  labs <- sort(unique(terra::values(ID2, mat = FALSE)))
  labs <- labs[is.finite(labs)]
  
  for (lab in labs) {
    m <- (ID2 == lab)
    p <- terra::patches(m, directions = 8)
    
    fr <- freq_noNA(p)
    if (is.null(fr) || nrow(fr) == 0) next
    
    lab_index <- as.integer(lab)
    keep_patch <- NA_real_
    
    if (lab_index >= 1 && lab_index <= length(pt_cells) && is.finite(pt_cells[lab_index])) {
      keep_patch <- terra::values(p, mat = FALSE)[pt_cells[lab_index]]
    }
    
    # Fallback: keep largest connected component
    if (!is.finite(keep_patch)) {
      fr <- fr[order(fr$count, decreasing = TRUE), , drop = FALSE]
      keep_patch <- fr$value[1]
    }
    
    # Drop ALL disconnected patches, not just small ones
    drop_patches <- fr$value[fr$value != keep_patch]
    
    if (length(drop_patches)) {
      ID2 <- terra::ifel(p %in% drop_patches, NA, ID2)
    }
  }
  
  # Fill removed cells only inside domain
  for (k in seq_len(max_iter_fill)) {
    vals <- terra::values(ID2, mat = FALSE)
    n_na <- sum(is.na(vals) & mask_vals == 1, na.rm = TRUE)
    if (n_na == 0) break
    
    neigh <- terra::focal(
      ID2,
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
    ID2 <- terra::setValues(ID2, vals)
  }
  
  # Final cleanup: remove any tiny stray fragments created during reassignment
  for (lab in labs) {
    m <- (ID2 == lab)
    p <- terra::patches(m, directions = 8)
    fr <- freq_noNA(p)
    if (is.null(fr) || nrow(fr) == 0) next
    
    lab_index <- as.integer(lab)
    keep_patch <- NA_real_
    
    if (lab_index >= 1 && lab_index <= length(pt_cells) && is.finite(pt_cells[lab_index])) {
      keep_patch <- terra::values(p, mat = FALSE)[pt_cells[lab_index]]
    }
    if (!is.finite(keep_patch)) {
      fr <- fr[order(fr$count, decreasing = TRUE), , drop = FALSE]
      keep_patch <- fr$value[1]
    }
    
    small_other <- fr$value[fr$value != keep_patch & fr$count < min_cells]
    if (length(small_other)) {
      ID2 <- terra::ifel(p %in% small_other, NA, ID2)
    }
  }
  
  vals <- terra::values(ID2, mat = FALSE)
  vals[mask_vals != 1] <- NA
  ID2 <- terra::setValues(ID2, vals)
  
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
#' @param boundary Optional `sf` polygon defining the tessellation domain. Used when `template_rast` is `NULL`.
#' @param template_rast Optional `terra::SpatRaster` template raster. Provide this instead of `boundary` + `res`.
#' @param method Character. Allocation method; one of `"argmin"` or `"partition"`.
#' @param weight_model Character. One of "multiplicative", "power", or "additive".
#'   Controls how distances and weights combine into effective cost.
#' @param weight_power Numeric > 0. Only used when weight_model = "power".
#'   Controls the distance exponent.


weighted_voronoi <- function(points_sf,
                             weight_col,
                             boundary = NULL,        # sf POLYGON/MULTIPOLYGON (optional)
                             template_rast = NULL,   # terra SpatRaster (optional)
                             res = NULL,             # numeric; only used if boundary is provided
                             weight_transform = function(w) w,
                             weight_model = c("multiplicative", "power", "additive"),
                             weight_power = 1,
                             method = c("argmin", "partition"),
                             max_dist = NULL,
                             verbose = TRUE,
                             island_min_cells = 5,
                             island_fill_iter = 50) {
  
  method <- match.arg(method)
  weight_model <- match.arg(weight_model)
  if (!is.finite(weight_power) || weight_power <= 0) {
    stop("weight_power must be finite and > 0.")
  }
  
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
    
    Ei <- .effective_cost(di, w[i], weight_model = weight_model, weight_power = weight_power)
    upd <- Ei < E
    E  <- terra::ifel(upd, Ei, E)
    ID <- terra::ifel(upd, i,  ID)
  }
  
  if (method == "argmin") {
    return(list(allocation = ID, cost_surface = E, weights = w, weights_raw = w_raw, weight_model = weight_model,
                weight_power = weight_power))
  }
  
  # ---- Partition-based output + island removal ----
  ID <- force_points_to_own_cell(ID, points_sf)
  
  domain_mask <- terra::ifel(!is.na(r), 1, NA)
  
  ID_clean <- remove_islands(
    ID,
    points_sf = points_sf,
    domain_mask = domain_mask,
    min_cells = island_min_cells,
    max_iter_fill = island_fill_iter
  )
  
  # Polygonise
  poly <- terra::as.polygons(ID_clean, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  poly_sf <- sf::st_as_sf(poly)
  names(poly_sf)[names(poly_sf) != "geometry"] <- "generator_id"
  
  cr <- terra::crs(r, proj = TRUE)
  if (!is.null(boundary_use)) {
    sf::st_crs(poly_sf) <- sf::st_crs(boundary_use)
  } else if (!is.na(cr) && nzchar(cr)) {
    sf::st_crs(poly_sf) <- sf::st_crs(cr)
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
    weights_raw = w_raw,
    weight_model = weight_model,
    weight_power = weight_power
  ))
}

#--------------------------------------------------------------------------------#

#' Weighted geodesic tessellation (core)
#'
#' Computes a weighted tessellation using domain-constrained (geodesic) distances.
#' Distances are calculated as shortest-path distances through a rasterised domain mask.
#'
#' @inheritParams weighted_voronoi_domain
#' @param weight_model Character. One of "multiplicative", "power", or "additive".
#'   Controls how distances and weights combine into effective cost.
#' @param weight_power Numeric > 0. Only used when weight_model = "power".
#'   Controls the distance exponent.
#' @param geodesic_engine Character. Geodesic allocation engine; one of
#'   `"classic"` or `"multisource"`.
#' @param return_polygons Logical. If `TRUE`, polygonise the cleaned allocation
#'   raster and attach point attributes. If `FALSE`, return allocation outputs only.
#' @param prepared Optional prepared geodesic context created by
#'   [prepare_geodesic_context()] for repeated compatible geodesic runs.
#' @return A list containing polygon output, allocation raster, and weights.
#' @export

weighted_voronoi_geodesic <- function(points_sf, weight_col, boundary_sf,
                                      res = 20,
                                      weight_transform = function(w) w,
                                      weight_model = c("multiplicative", "power", "additive"),
                                      weight_power = 1,
                                      close_mask = TRUE,
                                      close_iters = 1,
                                      resistance_rast = NULL,
                                      dem_rast = NULL,
                                      use_tobler = TRUE,
                                      tobler_v0_kmh = 6,
                                      tobler_a = 3.5,
                                      tobler_b = 0.05,
                                      min_speed_kmh = 0.25,
                                      anisotropy = c("none", "terrain"),
                                      uphill_factor = 1,
                                      downhill_factor = 1,
                                      island_min_cells = 5,
                                      island_fill_iter = 50,
                                      geodesic_engine = c("classic", "multisource"),
                                      return_polygons = TRUE,
                                      prepared = NULL,
                                      verbose = TRUE) {
  
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("raster", quietly = TRUE)) stop("Install raster.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  
  anisotropy <- match.arg(anisotropy)
  geodesic_engine <- match.arg(geodesic_engine)
  
  if (!is.null(prepared)) {
    if (!identical(sf::st_crs(boundary_sf), sf::st_crs(prepared$boundary_sf))) {
      stop("prepared context CRS does not match boundary_sf CRS.")
    }
  }
  
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
  
  # --- raster template + mask / prepared context ---
  if (!is.null(prepared)) {
    if (!is.list(prepared)) stop("prepared must be a prepared geodesic context object.")
    if (is.null(prepared$mask_r) || !inherits(prepared$mask_r, "SpatRaster")) {
      stop("prepared$mask_r must be a terra SpatRaster.")
    }
    if (is.null(prepared$transition)) {
      stop("prepared context does not contain a transition object.")
    }
    
    r <- prepared$mask_r
    tr <- prepared$transition
    
    prep_aniso <- prepared$settings$anisotropy
    prep_engine <- prepared$settings$geodesic_engine
    
    if (!identical(prep_aniso, anisotropy)) {
      stop("prepared context anisotropy does not match requested anisotropy.")
    }
    if (!identical(prep_engine, geodesic_engine)) {
      stop("prepared context geodesic_engine does not match requested geodesic_engine.")
    }
    
  } else {
    r <- .build_domain_mask(
      boundary_sf = boundary_sf,
      res = res,
      close_mask = close_mask,
      close_iters = close_iters
    )
  }
  
  # --- weights ---
  w_raw <- as.numeric(points_sf[[weight_col]])
  w <- weight_transform(w_raw)
  if (any(!is.finite(w)) || any(w <= 0)) stop("Weights must be finite and > 0.")
  
  weight_model <- match.arg(weight_model)
  if (!is.finite(weight_power) || weight_power <= 0) stop("weight_power must be finite and > 0.")
  
  # --- transition object ---
  if (anisotropy == "none") {
    resistance <- .prep_resistance(
      mask_r = r,
      resistance_rast = resistance_rast,
      dem_rast = dem_rast,
      use_tobler = use_tobler,
      tobler_v0_kmh = tobler_v0_kmh,
      tobler_a = tobler_a,
      tobler_b = tobler_b,
      min_speed_kmh = min_speed_kmh
    )
    
    tr <- .build_isotropic_transition(resistance)
    
  } else if (anisotropy == "terrain") {
    
    if (is.null(dem_rast)) {
      stop("anisotropy = 'terrain' requires dem_rast.")
    }
    
    tr <- .build_terrain_anisotropic_transition(
      dem_rast = dem_rast,
      mask_r = r,
      uphill_factor = uphill_factor,
      downhill_factor = downhill_factor,
      v0_kmh = tobler_v0_kmh,
      a = tobler_a,
      b = tobler_b,
      min_speed_kmh = min_speed_kmh
    )
  }
  
  if (geodesic_engine == "classic") {
    engine_out <- .compute_geodesic_allocation_classic(
      tr = tr,
      points_sf = points_sf,
      weight_vec = w,
      weight_model = weight_model,
      weight_power = weight_power,
      template_rast = r,
      verbose = verbose
    )
  } else if (geodesic_engine == "multisource") {
    engine_out <- .compute_geodesic_allocation_multisource(
      tr = tr,
      points_sf = points_sf,
      weight_vec = w,
      weight_model = weight_model,
      weight_power = weight_power,
      template_rast = r,
      graph = if (!is.null(prepared)) prepared$graph else NULL,
      anisotropy = anisotropy,
      verbose = verbose
    )
  }
  
  ID <- engine_out$allocation
  all_unreachable <- engine_out$unreachable
  
  domain_mask <- terra::ifel(!is.na(r), 1, NA)
  
  ID_clean <- .postprocess_allocation(
    ID = ID,
    points_sf = points_sf,
    domain_mask = domain_mask,
    island_min_cells = island_min_cells,
    island_fill_iter = island_fill_iter,
    fill_iter = 300
  )
  
  names(ID_clean) <- "allocation"
  
  poly_sf <- NULL
  
  if (return_polygons) {
    # Polygonise from cleaned allocation
    poly <- terra::as.polygons(ID_clean, dissolve = TRUE, values = TRUE, na.rm = TRUE)
    poly_sf <- sf::st_as_sf(poly)
    sf::st_crs(poly_sf) <- sf::st_crs(boundary_sf)
    names(poly_sf)[names(poly_sf) != "geometry"] <- "generator_id"
    
    # Final exact clipping to true boundary geometry
    poly_sf <- sf::st_make_valid(poly_sf)
    bnd_sf <- sf::st_make_valid(boundary_sf)
    poly_sf <- sf::st_intersection(poly_sf, bnd_sf)
    
    if (nrow(poly_sf) == 0) {
      stop("Geodesic clipping produced 0 polygons. Check boundary validity and CRS.")
    }
    
    poly_sf <- poly_sf |>
      dplyr::group_by(generator_id) |>
      dplyr::summarise(
        dplyr::across(dplyr::where(~ !inherits(.x, "sfc")), \(z) z[1]),
        do_union = TRUE
      ) |>
      dplyr::ungroup()
    
    # Attach attributes
    poly_sf$weight <- w[poly_sf$generator_id]
    pts_df <- sf::st_drop_geometry(points_sf)
    pts_df$generator_id <- seq_len(nrow(points_sf))
    poly_sf <- merge(poly_sf, pts_df, by = "generator_id", all.x = TRUE)
  }
  
  list(
    polygons = poly_sf,
    allocation = ID_clean,
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
#' @param dem_rast Optional SpatRaster providing elevation or resistance surface.
#'   Must align with the tessellation domain and resolution.
#' @param use_tobler Logical; if TRUE, apply Tobler's hiking function to convert
#'   slope into isotropic movement cost.
#' @param tobler_v0_kmh Base walking speed on flat terrain (km/h).
#' @param tobler_a Tobler exponential slope coefficient (default -3.5).
#' @param tobler_b Tobler slope multiplier (default 0.05).
#' @param min_speed_kmh Minimum allowed speed to avoid infinite costs.
#' @param verbose Logical. If `TRUE`, prints progress.
#' @param resistance_rast Optional SpatRaster giving movement resistance (>0). Overrides dem_rast/Tobler when provided.
#' @param weight_model Character. One of "multiplicative", "power", or "additive".
#'   Controls how distances and weights combine into effective cost.
#' @param weight_power Numeric > 0. Only used when weight_model = "power".
#'   Controls the distance exponent.
#' @param anisotropy Character. Directional cost model for geodesic distance.
#'   \describe{
#'     \item{"none"}{Standard isotropic geodesic distance (default).}
#'     \item{"terrain"}{Direction-dependent movement based on terrain slope (DEM required).}
#'   }
#'
#' @param uphill_factor Numeric > 0. Multiplier controlling additional cost of uphill movement
#'   when `anisotropy = "terrain"`. Values > 1 penalise uphill movement more strongly.
#'
#' @param downhill_factor Numeric > 0. Multiplier controlling ease of downhill movement
#'   when `anisotropy = "terrain"`. Values > 1 make downhill travel easier.
#'   
#' @param geodesic_engine Character. Geodesic allocation engine to use when
#'   `distance = "geodesic"`.
#'   \describe{
#'     \item{"classic"}{Per-generator accumulated-cost allocation. Supports all
#'     current geodesic modes and weight models.}
#'     \item{"multisource"}{Single-pass multisource allocation. Currently supported
#'     only for `weight_model = "additive"` and `anisotropy = "none"`.}
#'   } 
#' @param prepared Optional prepared geodesic context created by
#'   [prepare_geodesic_context()] for repeated compatible geodesic runs.
#'   
#' @details
#' When `distance = "geodesic"`, distances are computed as shortest paths
#' constrained to the spatial domain. If `dem_rast` is supplied and
#' `use_tobler = TRUE`, movement cost between adjacent raster cells is modified
#' using Tobler's hiking function, such that steeper slopes increase effective
#' distance. This allows elevation or resistance surfaces to influence spatial
#' allocation while preserving a complete tessellation.
#' When `distance = "geodesic"` and `anisotropy = "terrain"`, movement costs are
#' computed using a direction-dependent extension of a Tobler-like hiking function.
#'
#' Movement between raster cells becomes asymmetric: uphill and downhill transitions
#' have different costs. This results in anisotropic (direction-dependent) geodesic
#' tessellations.
#'
#' Currently, anisotropic terrain mode:
#' \itemize{
#'   \item requires a `dem_rast` input
#'   \item does not combine with a user-supplied `resistance_rast`
#'   \item uses 8-directional neighbourhood transitions
#' }
#' @return A list with elements including:
#' \describe{
#'   \item{polygons}{An `sf` object with one polygon per generator.}
#'   \item{allocation}{A `terra::SpatRaster` assigning each cell to a generator.}
#'   \item{summary}{A generator-level summary table.}
#'   \item{diagnostics}{A list of diagnostic metrics and settings.}
#' }
#' For geodesic allocation, `geodesic_engine = "classic"` computes one
#' accumulated-cost surface per generator and assigns each raster cell to the
#' minimum effective cost. This is the reference implementation and supports all
#' current geodesic modes.
#'
#' `geodesic_engine = "multisource"` provides a scalable alternative for
#' additive-weight isotropic geodesic tessellations. It uses a single multisource
#' shortest-path propagation and is currently available only when
#' `weight_model = "additive"` and `anisotropy = "none"`.
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
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' crs_use <- "EPSG:3857"
#'
#' boundary_sf <- st_sf(
#'   id = 1,
#'   geometry = st_sfc(
#'     st_polygon(list(rbind(
#'       c(0, 0), c(1000, 0), c(1000, 1000),
#'       c(0, 1000), c(0, 0)
#'     ))),
#'     crs = crs_use
#'   )
#' )
#'
#' points_sf <- st_sf(
#'   population = c(1, 1),
#'   geometry = st_sfc(
#'     st_point(c(200, 500)),
#'     st_point(c(800, 500)),
#'     crs = crs_use
#'   )
#' )
#'
#' dem_rast <- rast(
#'   ext = ext(0, 1000, 0, 1000),
#'   resolution = 50,
#'   crs = crs_use
#' )
#'
#' xy <- crds(dem_rast, df = TRUE)
#' values(dem_rast) <- xy$x * 20
#'
#' out <- weighted_voronoi_domain(
#'   points_sf = points_sf,
#'   weight_col = "population",
#'   boundary_sf = boundary_sf,
#'   distance = "geodesic",
#'   dem_rast = dem_rast,
#'   anisotropy = "terrain",
#'   uphill_factor = 3,
#'   downhill_factor = 1.2
#' )
#' }
#' @export

weighted_voronoi_domain <- function(points_sf,
                                    weight_col,
                                    boundary_sf,
                                    res = 20,
                                    weight_transform = function(w) w,
                                    weight_model = c("multiplicative", "power", "additive"),
                                    weight_power = 1,
                                    distance = c("euclidean", "geodesic"),
                                    # Euclidean options
                                    max_dist = NULL,
                                    island_min_cells = 5,
                                    island_fill_iter = 50,
                                    clip_to_boundary = TRUE,
                                    # Geodesic options
                                    close_mask = TRUE,
                                    close_iters = 1,
                                    resistance_rast = NULL,
                                    dem_rast = NULL,
                                    use_tobler = TRUE,
                                    tobler_v0_kmh = 6,
                                    tobler_a = 3.5,
                                    tobler_b = 0.05,
                                    min_speed_kmh = 0.25,
                                    anisotropy = c("none", "terrain"),
                                    uphill_factor = 1,
                                    downhill_factor = 1,
                                    geodesic_engine = c("classic", "multisource"),
                                    prepared = NULL,
                                    # general
                                    verbose = TRUE) {
  
  distance <- match.arg(distance)
  anisotropy <- match.arg(anisotropy)
  geodesic_engine <- match.arg(geodesic_engine)
  
  if (!is.null(prepared) && distance != "geodesic") {
    stop("prepared can only be used when distance = 'geodesic'.")
  }
  
  if (!anisotropy %in% c("none", "terrain")) {
    stop("anisotropy must be 'none' or 'terrain'.")
  }
  
  if (anisotropy == "terrain") {
    if (is.null(dem_rast)) {
      stop("anisotropy = 'terrain' requires dem_rast.")
    }
    if (!is.null(resistance_rast)) {
      warning("resistance_rast is currently ignored when anisotropy = 'terrain'.")
    }
  }
  
  if (!is.finite(uphill_factor) || uphill_factor <= 0) {
    stop("uphill_factor must be finite and > 0.")
  }
  
  if (!is.finite(downhill_factor) || downhill_factor <= 0) {
    stop("downhill_factor must be finite and > 0.")
  }
  
  if (geodesic_engine == "multisource" && weight_model != "additive") {
    stop("geodesic_engine = 'multisource' currently requires weight_model = 'additive'.")
  }
  
  if (geodesic_engine == "multisource" && anisotropy != "none") {
    stop("geodesic_engine = 'multisource' currently requires anisotropy = 'none'.")
  }
  
  weight_model <- match.arg(weight_model)
  if (!is.finite(weight_power) || weight_power <= 0) stop("weight_power must be finite and > 0.")
  
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
      weight_model = weight_model,
      weight_power = weight_power,
      method = "partition",
      max_dist = max_dist,
      verbose = verbose,
      island_min_cells = island_min_cells,
      island_fill_iter = island_fill_iter
    )
    
    if (clip_to_boundary && !is.null(boundary_sf)) {
      poly_sf <- sf::st_make_valid(out$polygons)
      b <- sf::st_make_valid(boundary_sf)
      
      poly_sf <- sf::st_intersection(poly_sf, b)
      if (nrow(poly_sf) == 0) stop("Clipping produced 0 polygons. Check boundary CRS/validity.")
      
      # Dissolve back to one feature per generator_id, KEEPING attributes
      poly_sf <- poly_sf |>
        dplyr::group_by(generator_id) |>
        dplyr::summarise(
          dplyr::across(dplyr::where(~ !inherits(.x, "sfc")), \(z) z[1]),
          do_union = TRUE
        ) |>
        dplyr::ungroup()
      
      out$polygons <- poly_sf
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
        clip_to_boundary = clip_to_boundary,
        weight_model = weight_model,
        weight_power = weight_power
      )
    )
    
    return(list(
      method = "euclidean",
      weight_model = weight_model,
      weight_power = weight_power,
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
    weight_model = weight_model,
    weight_power = weight_power,
    close_mask = close_mask,
    close_iters = close_iters,
    dem_rast = dem_rast,
    resistance_rast = resistance_rast,
    use_tobler = use_tobler,
    tobler_v0_kmh = tobler_v0_kmh,
    tobler_a = tobler_a,
    tobler_b = tobler_b,
    min_speed_kmh = min_speed_kmh,
    anisotropy = anisotropy,
    uphill_factor = uphill_factor,
    downhill_factor = downhill_factor,
    island_min_cells = island_min_cells,
    island_fill_iter = island_fill_iter,
    geodesic_engine = geodesic_engine,
    prepared = prepared,
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
      close_iters = close_iters,
      dem_rast = dem_rast,
      use_tobler = use_tobler,
      tobler_v0_kmh = tobler_v0_kmh,
      tobler_a = tobler_a,
      tobler_b = tobler_b,
      min_speed_kmh = min_speed_kmh,
      weight_model = weight_model,
      weight_power = weight_power,
      island_min_cells = island_min_cells,
      island_fill_iter = island_fill_iter,
      resistance_rast_supplied = !is.null(resistance_rast),
      dem_rast_supplied = !is.null(dem_rast),
      anisotropy = anisotropy,
      uphill_factor = uphill_factor,
      downhill_factor = downhill_factor,
      geodesic_engine = geodesic_engine,
      prepared_supplied = !is.null(prepared)
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
    diagnostics = diagnostics,
    weight_model = weight_model,
    weight_power = weight_power
  )
}

