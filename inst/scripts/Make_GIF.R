#install.packages(c("terra", "sf", "magick"))

library(sf)

crs_use <- 32636  # any projected CRS in metres

boundary_rect <- st_sf(
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0),
    c(1200, 0),
    c(1200, 800),
    c(0, 800),
    c(0, 0)
  )))),
  crs = crs_use
)

set.seed(1)
pts <- st_sample(boundary_rect, size = 8)
points_rect <- st_sf(
  id = 1:length(pts),
  population = round(exp(rnorm(length(pts), log(300), 0.7))),
  geometry = pts
)

library(terra)
library(sf)
library(magick)

make_euclidean_tessellation_gif <- function(points_sf,
                                            weight_col,
                                            boundary_sf,
                                            res = 10,
                                            weight_transform = log10,
                                            n_frames = 25,
                                            out_gif = "man/figures/animation-euclidean.gif",
                                            width_px = 900,
                                            height_px = 600,
                                            dpi = 120) {
  stopifnot(inherits(points_sf, "sf"), inherits(boundary_sf, "sf"))
  if (!weight_col %in% names(points_sf)) stop("weight_col not found in points_sf")
  
  # CRS harmonise
  if (sf::st_crs(points_sf) != sf::st_crs(boundary_sf)) {
    boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(points_sf))
  }
  if (sf::st_is_longlat(points_sf)) stop("Use a projected CRS in metres.")
  
  boundary_sf <- sf::st_make_valid(boundary_sf)
  
  # Raster domain mask
  bnd_v <- terra::vect(boundary_sf)
  r <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
  r <- terra::setValues(r, 1)
  r <- terra::mask(r, bnd_v)
  
  # Weights
  w_raw <- as.numeric(points_sf[[weight_col]])
  w <- weight_transform(w_raw)
  if (any(!is.finite(w)) || any(w <= 0)) stop("Weights must be finite and > 0 after transform.")
  
  pts_v <- terra::vect(points_sf)
  
  # Build weighted cost stack: (distance / weight)
  cost_list <- vector("list", nrow(points_sf))
  for (i in seq_len(nrow(points_sf))) {
    di <- terra::distance(r, pts_v[i])
    cost_list[[i]] <- di / w[i]
  }
  S <- terra::rast(cost_list)
  
  # Final allocation (full tessellation)
  ID_final <- terra::which.min(S)
  
  # Minimum cost surface (for revealing animation)
  min_cost <- terra::app(S, fun = function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(NA_real_)
    min(vals)
  })
  
  # --- FIXED: compute thresholds from values (base quantile) ---
  q <- seq(0.02, 1, length.out = n_frames)
  vmin <- terra::values(min_cost, mat = FALSE)
  vmin <- vmin[is.finite(vmin)]
  if (!length(vmin)) stop("min_cost has no finite values inside the domain (check boundary/raster mask).")
  th <- as.numeric(stats::quantile(vmin, probs = q, na.rm = TRUE, names = FALSE))
  
  # Ensure output folder exists
  dir.create(dirname(out_gif), recursive = TRUE, showWarnings = FALSE)
  
  # Make frames
  frame_files <- character(n_frames)
  for (k in seq_len(n_frames)) {
    # show only cells whose min_cost <= threshold; others are NA (blank)
    ID_k <- terra::ifel(min_cost <= th[k], ID_final, NA)
    
    f <- file.path(tempdir(), sprintf("frame_%03d.png", k))
    frame_files[k] <- f
    
    png(f, width = width_px, height = height_px, res = dpi)
    
    plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2,
         main = "Weighted Euclidean tessellation (revealed)")
    plot(ID_k, add = TRUE, axes = FALSE, legend = FALSE)
    plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red", col = "black", cex = 1.1)
    text(sf::st_coordinates(points_sf), labels = points_sf[[weight_col]], pos = 3, cex = 0.8)
    
    dev.off()
  }
  
  # Assemble GIF
  img <- magick::image_read(frame_files)
  img <- magick::image_animate(img, fps = 10)
  magick::image_write(img, path = out_gif)
  
  message("Wrote GIF: ", out_gif)
  invisible(out_gif)
}


make_euclidean_tessellation_gif(
  points_sf = points_rect,
  weight_col = "population",
  boundary_sf = boundary_rect,
  res = 2,
  weight_transform = log10,
  n_frames = 25,
  out_gif = "man/figures/animation-euclidean-four.gif"
)
