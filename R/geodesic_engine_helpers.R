.compute_geodesic_allocation_classic <- function(
    tr,
    points_sf,
    weight_vec,
    weight_model,
    weight_power,
    template_rast = NULL,
    verbose = TRUE
) {
  S <- .compute_geodesic_cost_stack(
    tr = tr,
    points_sf = points_sf,
    weight_vec = weight_vec,
    weight_model = weight_model,
    weight_power = weight_power,
    verbose = verbose
  )
  
  all_unreachable <- .detect_all_unreachable(S)
  
  ID <- .assign_from_cost_stack(
    cost_stack = S,
    all_unreachable = all_unreachable
  )
  
  list(
    allocation = ID,
    unreachable = all_unreachable,
    cost_stack = S
  )
}

.compute_geodesic_allocation_multisource <- function(
    tr,
    points_sf,
    weight_vec,
    weight_model,
    weight_power,
    template_rast,
    anisotropy = "none",
    verbose = TRUE
) {
  if (!requireNamespace("gdistance", quietly = TRUE)) stop("Install gdistance.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install terra.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix.")
  
  if (weight_model != "additive") {
    stop("geodesic_engine = 'multisource' currently requires weight_model = 'additive'.")
  }
  if (!identical(anisotropy, "none")) {
    stop("geodesic_engine = 'multisource' currently requires anisotropy = 'none'.")
  }
  if (is.null(template_rast) || !inherits(template_rast, "SpatRaster")) {
    stop("template_rast must be supplied as a terra SpatRaster.")
  }
  
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
  
  # Source cells from point locations
  pts_v <- terra::vect(points_sf)
  xy <- terra::crds(pts_v)
  src_cells <- terra::cellFromXY(template_rast, xy)
  
  if (any(!is.finite(src_cells))) {
    stop("Could not map one or more points to raster cells.")
  }
  
  # Distances and owner labels
  dist <- rep(Inf, n)
  owner <- rep(NA_integer_, n)
  
  # Valid domain cells only
  mask_vals <- terra::values(template_rast, mat = FALSE)
  inside_domain <- !is.na(mask_vals)
  
  # Initialise each source with additive offset = 1 / weight
  # If multiple generators fall in the same cell, keep the lower initial cost;
  # break exact ties by smaller generator index for determinism.
  for (i in seq_along(src_cells)) {
    cell <- src_cells[i]
    d0 <- 1 / weight_vec[i]
    
    if (!inside_domain[cell]) next
    
    if (d0 < dist[cell] || (isTRUE(all.equal(d0, dist[cell])) && (is.na(owner[cell]) || i < owner[cell]))) {
      dist[cell] <- d0
      owner[cell] <- i
    }
  }
  
  # Priority queue (min-heap) with lazy deletion
  heap <- .heap_new(max(1024L, length(src_cells) * 4L))
  
  seeded <- which(is.finite(dist) & !is.na(owner))
  for (cell in seeded) {
    .heap_push(heap, key = dist[cell], node = cell)
  }
  
  tol <- 1e-12
  
  n_pop <- 0L
  while (heap$size > 0L) {
    item <- .heap_pop(heap)
    u <- item$node
    du <- item$key
    
    # Skip stale heap entries
    if (!is.finite(dist[u]) || abs(du - dist[u]) > tol) next
    
    if (verbose) {
      n_pop <- n_pop + 1L
      if (n_pop %% 50000L == 0L) {
        message(sprintf("Multisource Dijkstra settled %d nodes", n_pop))
      }
    }
    
    idx1 <- start_idx[u]
    idx2 <- start_idx[u + 1L] - 1L
    if (idx2 < idx1) next
    
    nbrs <- idx1:idx2
    
    for (k in nbrs) {
      v <- to[k]
      nd <- du + edge_cost[k]
      cand_owner <- owner[u]
      
      # Standard relaxation
      if (nd + tol < dist[v]) {
        dist[v] <- nd
        owner[v] <- cand_owner
        .heap_push(heap, key = nd, node = v)
        
      } else if (abs(nd - dist[v]) <= tol) {
        # Deterministic tie-break: prefer smaller owner id
        if (!is.na(cand_owner) && (is.na(owner[v]) || cand_owner < owner[v])) {
          owner[v] <- cand_owner
          .heap_push(heap, key = nd, node = v)
        }
      }
    }
  }
  
  # Outside-domain cells should remain NA
  owner[!inside_domain] <- NA_integer_
  
  allocation <- template_rast
  allocation <- terra::setValues(allocation, owner)
  
  unreachable <- template_rast
  unreachable_vals <- inside_domain & !is.finite(dist)
  unreachable_vals[!inside_domain] <- NA
  unreachable <- terra::setValues(unreachable, unreachable_vals)
  
  list(
    allocation = allocation,
    unreachable = unreachable,
    distances = dist
  )
}

.heap_new <- function(capacity = 1024L) {
  e <- new.env(parent = emptyenv())
  e$keys <- rep(NA_real_, capacity)
  e$nodes <- rep(NA_integer_, capacity)
  e$size <- 0L
  e
}

.heap_push <- function(heap, key, node) {
  if (heap$size >= length(heap$keys)) {
    new_cap <- max(2L * length(heap$keys), 1L)
    length(heap$keys) <- new_cap
    length(heap$nodes) <- new_cap
  }
  
  i <- heap$size + 1L
  heap$size <- i
  
  while (i > 1L) {
    parent <- i %/% 2L
    if (!is.na(heap$keys[parent]) && heap$keys[parent] <= key) break
    heap$keys[i] <- heap$keys[parent]
    heap$nodes[i] <- heap$nodes[parent]
    i <- parent
  }
  
  heap$keys[i] <- key
  heap$nodes[i] <- node
  invisible(NULL)
}

.heap_pop <- function(heap) {
  if (heap$size < 1L) stop("Heap is empty.")
  
  out_key <- heap$keys[1L]
  out_node <- heap$nodes[1L]
  
  last_key <- heap$keys[heap$size]
  last_node <- heap$nodes[heap$size]
  heap$size <- heap$size - 1L
  
  if (heap$size >= 1L) {
    i <- 1L
    repeat {
      left <- 2L * i
      right <- left + 1L
      if (left > heap$size) break
      
      child <- left
      if (right <= heap$size && heap$keys[right] < heap$keys[left]) {
        child <- right
      }
      
      if (last_key <= heap$keys[child]) break
      
      heap$keys[i] <- heap$keys[child]
      heap$nodes[i] <- heap$nodes[child]
      i <- child
    }
    
    heap$keys[i] <- last_key
    heap$nodes[i] <- last_node
  }
  
  list(key = out_key, node = out_node)
}