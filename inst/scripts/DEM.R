
library(sf)
library(terra)
library(weightedVoronoi)

# 1) Use a projected CRS (metres)
crs_use <- 32636

# 2) Define a rectangular boundary (fast + robust for demos)
boundary_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0),
    c(1200, 0),
    c(1200, 800),
    c(0, 800),
    c(0, 0)
  )))),
  crs = crs_use
)

# 3) Sample generator points INSIDE boundary (guaranteed inside)
set.seed(1)
n_pts <- 8
pts <- st_sample(boundary_sf, size = n_pts, type = "random")

points_sf <- st_sf(
  id = seq_len(length(pts)),
  population = round(exp(rnorm(length(pts), log(300), 0.8))),  # skewed weights
  geometry = pts,
  crs = crs_use
)

# Sanity check (should be all TRUE)
inside <- st_within(points_sf, boundary_sf, sparse = FALSE)[, 1]
stopifnot(all(inside))

# 4) Choose resolution (cellsize) — keep consistent everywhere
cellsize <- 10
bnd_v <- terra::vect(boundary_sf)

# 5) Create a synthetic DEM with a strong ridge (easy to see)
dem <- terra::rast(
  ext = terra::ext(bnd_v),
  resolution = cellsize,
  crs = terra::crs(bnd_v)
)

xy <- terra::crds(dem, df = TRUE)

# Ridge parameters (tweak height to strengthen/weaken effect)
y0 <- (min(xy$y) + max(xy$y)) / 2
sigma <- (max(xy$y) - min(xy$y)) / 14
height <- 1000

terra::values(dem) <- height * exp(-((xy$y - y0)^2) / (2 * sigma^2))

par(mfrow = c(1,1))
# Optional plot check
plot(st_geometry(boundary_sf), border = "black", lwd = 2, main = "Demo inputs")
plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")


out_geo_plain <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  verbose = FALSE
)

out_geo_dem <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  dem_rast = dem,
  use_tobler = TRUE,
  verbose = FALSE
)


plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_geo_plain$polygons["generator_id"], main = "Geodesic", add = TRUE)
plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")

plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2)
plot(out_geo_dem$polygons["generator_id"], main = "Geodesic + Tobler ridge resistance", add = TRUE)
plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red")


plot(sf::st_geometry(boundary_sf), border = "black", lwd = 2,
     main = "Tessellation over DEM contours")

# Contours from DEM
terra::contour(dem, add = TRUE, nlevels = 6, col = "grey50", lwd = 0.7)

# Polygons + points
plot(out_geo_dem$polygons["generator_id"], add = TRUE, border = "white", lwd = 0.4, reset = FALSE)
plot(sf::st_geometry(points_sf), add = TRUE, pch = 21, bg = "red", col = "black")

library(sf)
library(terra)
library(ggplot2)
library(png)
library(rayshader)

# ---------- helpers ----------

dem_to_matrix <- function(dem_rast) {
  stopifnot(inherits(dem_rast, "SpatRaster"))
  m <- terra::as.matrix(dem_rast, wide = TRUE)
  # flip vertically to match rayshader orientation
  m[nrow(m):1, ]
}

save_tess_overlay_png <- function(polygons_sf, boundary_sf, points_sf,
                                  file,
                                  width = 1600, height = 1100, dpi = 220,
                                  alpha_polys = 0.65) {
  stopifnot(inherits(polygons_sf, "sf"), inherits(boundary_sf, "sf"), inherits(points_sf, "sf"))
  
  # enforce same bbox so overlay aligns with DEM exactly
  bb <- sf::st_bbox(boundary_sf)
  
  p <- ggplot() +
    geom_sf(data = polygons_sf, aes(fill = factor(generator_id)),
            colour = "white", linewidth = 0.2, alpha = alpha_polys) +
    geom_sf(data = boundary_sf, fill = NA, colour = "black", linewidth = 0.7) +
    geom_sf(data = points_sf, shape = 21, fill = "red", colour = "black",
            size = 2.3, stroke = 0.35) +
    guides(fill = "none") +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_void()
  
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, p, width = width / dpi, height = height / dpi, dpi = dpi, bg = "transparent")
  file
}


#----------------------------------------------------------------------------#

make_demo_dem <- function(boundary_sf, cellsize = 5, seed = 1,
                          n_hills = 9, hill_height = 1900,
                          ridge_height = 400, ridge_sigma_frac = 18,
                          smooth = 4) {
  set.seed(seed)
  bnd_v <- terra::vect(boundary_sf)
  
  dem <- terra::rast(ext = terra::ext(bnd_v), resolution = cellsize, crs = terra::crs(bnd_v))
  xy  <- terra::crds(dem, df = TRUE)
  
  # base: gentle slope (helps shading)
  z <- 0.08 * (xy$x - min(xy$x)) + 0.03 * (xy$y - min(xy$y))
  
  # add several broad hills
  for (k in seq_len(n_hills)) {
    cx <- runif(1, min(xy$x), max(xy$x))
    cy <- runif(1, min(xy$y), max(xy$y))
    sig <- runif(1,
                 (max(xy$x) - min(xy$x)) / 10,
                 (max(xy$x) - min(xy$x)) / 5)
    amp <- runif(1, 0.4, 1.0) * hill_height
    z <- z + amp * exp(-((xy$x - cx)^2 + (xy$y - cy)^2) / (2 * sig^2))
  }
  
  # add a ridge (N–S-ish band)
  y0 <- (min(xy$y) + max(xy$y)) / 2
  sigma <- (max(xy$y) - min(xy$y)) / ridge_sigma_frac
  z <- z + ridge_height * exp(-((xy$y - y0)^2) / (2 * sigma^2))
  
  terra::values(dem) <- z
  
  # mask to boundary & smooth a bit (makes it less pixel-y)
  dem <- terra::mask(dem, bnd_v)
  
  if (smooth > 0) {
    w <- matrix(1, 2 * smooth + 1, 2 * smooth + 1)
    dem <- terra::focal(dem, w = w, fun = mean, na.policy = "omit", fillvalue = NA)
  }
  
  dem
}


cellsize <- 2  # higher-res DEM (was 10)
dem <- make_demo_dem(boundary_sf, cellsize = cellsize, seed = 1, smooth = 3)
plot(dem); plot(terra::vect(boundary_sf), add = TRUE)


cellsize <- 2  # try 3 or 5 first; 1 can get heavy fast
out_geo_dem <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  dem_rast = dem,
  use_tobler = TRUE,
  verbose = FALSE
)

overlay_file <- save_tess_overlay_png(
  polygons_sf = out_geo_dem$polygons,
  boundary_sf = boundary_sf,
  points_sf   = points_sf,
  file        = "man/figures/fig-dem-overlay.png",
  width = 2600, height = 1700, dpi = 300,
  alpha_polys = 0.6
)

elmat <- dem_to_matrix(dem)
overlay_img <- png::readPNG("man/figures/fig-dem-overlay.png")

elmat %>%
  rayshader::height_shade() %>%
  rayshader::add_overlay(overlay_img, alphalayer = 1) %>%
  rayshader::plot_3d(
    elmat,
    zscale = 40,      # try 25–60 depending on DEM amplitude
    solid = TRUE,
    shadow = TRUE,
    phi = 30, theta = 45, zoom = 0.6,
    windowsize = c(1400, 1000)
  )

rayshader::render_snapshot(
  "man/figures/fig-dem-3d.png",
  title_text = "Geodesic tessellation with elevation resistance",
  title_size = 24
)
rgl::rgl.close()


#-------------------------------------------------------------------------------#

out_geo_flat <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "geodesic",
  # crucial: no dem_rast, no tobler
  verbose = FALSE
)

# flat surface on same grid as your DEM
bnd_v <- terra::vect(boundary_sf)

dem_flat <- terra::rast(
  ext = terra::ext(bnd_v),
  resolution = cellsize,
  crs = terra::crs(bnd_v)
)
dem_flat <- terra::setValues(dem_flat, 0)
dem_flat <- terra::mask(dem_flat, bnd_v)

elmat_flat <- dem_to_matrix(dem_flat)


flat_overlay_file <- save_tess_overlay_png(
  polygons_sf = out_geo_flat$polygons,
  boundary_sf = boundary_sf,
  points_sf   = points_sf,
  file        = "man/figures/fig-flat-overlay.png",
  width = 2600, height = 1700, dpi = 300,
  alpha_polys = 0.6
)


library(png)
library(rayshader)

flat_overlay_img <- png::readPNG("man/figures/fig-flat-overlay.png")

# base shading on flat surface
hs_flat <- rayshader::height_shade(elmat_flat)

png("man/figures/fig-flat-2d.png", width = 2600, height = 1700, res = 300)
hs_flat %>%
  rayshader::add_overlay(flat_overlay_img, alphalayer = 1) %>%
  rayshader::plot_map()
dev.off()

message("Wrote: man/figures/fig-flat-2d.png")


elmat_flat %>%
  rayshader::height_shade() %>%
  rayshader::add_overlay(flat_overlay_img, alphalayer = 1) %>%
  rayshader::plot_3d(
    elmat_flat,
    zscale = 40,         # doesn’t matter much because it’s flat
    solid = TRUE,
    shadow = TRUE,
    phi = 30, theta = 45, zoom = 0.60,
    windowsize = c(1400, 1000)
  )

rayshader::render_snapshot("man/figures/fig-flat-3d.png")
rgl::rgl.close()

message("Wrote: man/figures/fig-flat-3d.png")

#-------------------------------------------------------------------------------#

cellsize <- 2  # try 3 or 5 first; 1 can get heavy fast
out_geo_Euc <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = cellsize,
  weight_transform = log10,
  distance = "euclidean",
  verbose = FALSE
)

overlay_file <- save_tess_overlay_png(
  polygons_sf = out_geo_Euc$polygons,
  boundary_sf = boundary_sf,
  points_sf   = points_sf,
  file        = "man/figures/fig-Euc-overlay.png",
  width = 2600, height = 1700, dpi = 300,
  alpha_polys = 0.6
)

elmat <- dem_to_matrix(dem)
overlay_img <- png::readPNG("man/figures/fig-dem-overlay.png")