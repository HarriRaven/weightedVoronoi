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
    template_rast = template_rast,
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
    tr = NULL,
    points_sf,
    weight_vec,
    weight_model,
    weight_power,
    template_rast = NULL,
    graph = NULL,
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
  if (is.null(graph)) {
    if (is.null(template_rast) || !inherits(template_rast, "SpatRaster")) {
      stop("template_rast must be supplied as a terra SpatRaster when graph is NULL.")
    }
    if (is.null(tr)) {
      stop("tr must be supplied when graph is NULL.")
    }
  }
  
  if (!is.null(graph)) {
    from <- graph$from
    to <- graph$to
    edge_cost <- graph$edge_cost
    start_idx <- graph$start_idx
    template_rast <- graph$template_rast
  } else {
    graph <- .prepare_multisource_graph(tr = tr, template_rast = template_rast)
    from <- graph$from
    to <- graph$to
    edge_cost <- graph$edge_cost
    start_idx <- graph$start_idx
    template_rast <- graph$template_rast
  }
  
  n <- terra::ncell(template_rast)
  
  # --- DEBUG BLOCK 1: immediately after graph/template setup ---
  mask_vals <- terra::values(template_rast, mat = FALSE)
  inside_mask <- !is.na(mask_vals)
  
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
  
  # Cells actually represented in the multisource graph
  inside_graph <- rep(FALSE, n)
  graph_nodes <- unique(c(from, to))
  inside_graph[graph_nodes] <- TRUE
  
  # Effective domain for multisource propagation
  inside_domain <- inside_mask & inside_graph
  
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
  
  alloc_vals <- terra::values(allocation, mat = FALSE)
  
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