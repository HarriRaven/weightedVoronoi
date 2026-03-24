.prepare_multisource_graph <- function(tr, template_rast) {
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix.")
  if (!inherits(template_rast, "SpatRaster")) stop("template_rast must be a terra SpatRaster.")

  # Transition matrix stores conductance; convert to edge cost
  TM <- gdistance::transitionMatrix(tr)
  sm <- Matrix::summary(TM)

  # Keep only finite positive conductance edges
  ok <- is.finite(sm$x) & sm$x > 0
  if (!any(ok)) stop("Transition matrix has no finite positive edges.")

  from0 <- sm$i[ok]
  to0   <- sm$j[ok]
  cond0 <- sm$x[ok]

  # For isotropic multisource mode, always make the graph explicitly bidirectional.
  # This avoids any ambiguity about how the sparse matrix stores symmetric entries.
  offdiag <- from0 != to0

  from <- c(from0, to0[offdiag])
  to   <- c(to0,   from0[offdiag])
  conductance <- c(cond0, cond0[offdiag])

  edge_cost <- 1 / conductance

  # Sort by source node to allow fast outgoing-edge lookup
  ord <- order(from, to)
  from <- from[ord]
  to <- to[ord]
  edge_cost <- edge_cost[ord]

  n <- terra::ncell(template_rast)
  counts <- tabulate(from, nbins = n)
  start_idx <- cumsum(c(1L, counts))  # outgoing edges for u are start_idx[u]:(start_idx[u+1]-1)

  list(
    from = from,
    to = to,
    edge_cost = edge_cost,
    start_idx = start_idx,
    template_rast = template_rast
  )
}

.prepare_geodesic_context <- function(
    boundary_sf,
    res = 20,
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
    geodesic_engine = c("classic", "multisource")
) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf.")
  if (!inherits(boundary_sf, "sf")) stop("boundary_sf must be an sf object.")
  
  anisotropy <- match.arg(anisotropy)
  geodesic_engine <- match.arg(geodesic_engine)
  
  boundary_sf <- sf::st_make_valid(boundary_sf)
  
  mask_r <- .build_domain_mask(
    boundary_sf = boundary_sf,
    res = res,
    close_mask = close_mask,
    close_iters = close_iters
  )
  
  tr <- NULL
  graph <- NULL
  resistance <- NULL
  
  if (anisotropy == "none") {
    resistance <- .prep_resistance(
      mask_r = mask_r,
      resistance_rast = resistance_rast,
      dem_rast = dem_rast,
      use_tobler = use_tobler,
      tobler_v0_kmh = tobler_v0_kmh,
      tobler_a = tobler_a,
      tobler_b = tobler_b,
      min_speed_kmh = min_speed_kmh
    )
    
    tr <- .build_isotropic_transition(resistance)
    
    if (geodesic_engine == "multisource") {
      graph <- .prepare_multisource_graph(tr = tr, template_rast = mask_r)
    }
    
  } else {
    if (is.null(dem_rast)) {
      stop("anisotropy = 'terrain' requires dem_rast.")
    }
    
    tr <- .build_terrain_anisotropic_transition(
      dem_rast = dem_rast,
      mask_r = mask_r,
      uphill_factor = uphill_factor,
      downhill_factor = downhill_factor,
      v0_kmh = tobler_v0_kmh,
      a = tobler_a,
      b = tobler_b,
      min_speed_kmh = min_speed_kmh
    )
  }
  
  list(
    boundary_sf = boundary_sf,
    mask_r = mask_r,
    resistance = resistance,
    transition = tr,
    graph = graph,
    settings = list(
      res = res,
      close_mask = close_mask,
      close_iters = close_iters,
      anisotropy = anisotropy,
      geodesic_engine = geodesic_engine,
      use_tobler = use_tobler,
      tobler_v0_kmh = tobler_v0_kmh,
      tobler_a = tobler_a,
      tobler_b = tobler_b,
      min_speed_kmh = min_speed_kmh,
      uphill_factor = uphill_factor,
      downhill_factor = downhill_factor,
      resistance_rast_supplied = !is.null(resistance_rast),
      dem_rast_supplied = !is.null(dem_rast)
    )
  )
}

#' Prepare a geodesic context for repeated runs
#'
#' Precomputes the domain mask, aligned resistance handling, transition object,
#' and multisource graph representation (when applicable) for repeated geodesic
#' tessellation workflows.
#'
#' @param boundary_sf An `sf` POLYGON/MULTIPOLYGON defining the domain.
#' @param res Numeric. Raster resolution in CRS units (e.g. metres).
#' @param close_mask Logical. If `TRUE`, applies a morphological closing to the raster mask.
#' @param close_iters Integer. Number of closing iterations.
#' @param resistance_rast Optional SpatRaster giving movement resistance (>0).
#' @param dem_rast Optional SpatRaster providing elevation or resistance surface.
#' @param use_tobler Logical; if `TRUE`, apply Tobler's hiking function to convert slope into isotropic movement cost.
#' @param tobler_v0_kmh Base walking speed on flat terrain (km/h).
#' @param tobler_a Tobler exponential slope coefficient.
#' @param tobler_b Tobler slope multiplier.
#' @param min_speed_kmh Minimum allowed speed to avoid infinite costs.
#' @param anisotropy Character. One of `"none"` or `"terrain"`.
#' @param uphill_factor Numeric > 0. Additional uphill movement penalty when `anisotropy = "terrain"`.
#' @param downhill_factor Numeric > 0. Relative ease of downhill movement when `anisotropy = "terrain"`.
#' @param geodesic_engine Character. One of `"classic"` or `"multisource"`.
#'
#' @return A prepared geodesic context object for repeated geodesic allocation.
#' @export
prepare_geodesic_context <- function(
    boundary_sf,
    res = 20,
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
    geodesic_engine = c("classic", "multisource")
) {
  .prepare_geodesic_context(
    boundary_sf = boundary_sf,
    res = res,
    close_mask = close_mask,
    close_iters = close_iters,
    resistance_rast = resistance_rast,
    dem_rast = dem_rast,
    use_tobler = use_tobler,
    tobler_v0_kmh = tobler_v0_kmh,
    tobler_a = tobler_a,
    tobler_b = tobler_b,
    min_speed_kmh = min_speed_kmh,
    anisotropy = anisotropy,
    uphill_factor = uphill_factor,
    downhill_factor = downhill_factor,
    geodesic_engine = geodesic_engine
  )
}