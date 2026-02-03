#' Tobler isotropic resistance from DEM
#'
#' Builds a resistance raster (seconds per metre) from a DEM using Tobler's hiking function.
#' Isotropic: uses slope magnitude (no uphill vs downhill asymmetry).
#'
#' @param dem_rast terra SpatRaster DEM (elevation, metres)
#' @param mask_r terra SpatRaster mask (1 inside domain, NA outside)
#' @param v0_kmh base speed (km/h), default 6
#' @param a Tobler steepness coefficient, default 3.5
#' @param b Tobler offset, default 0.05
#' @param min_speed_kmh minimum speed clamp (km/h) to avoid infinite resistance
#' @return terra SpatRaster resistance (s/m), masked to domain
#' @noRd
.tobler_resistance_from_dem <- function(dem_rast, mask_r,
                                        v0_kmh = 6,
                                        a = 3.5,
                                        b = 0.05,
                                        min_speed_kmh = 0.25) {
  stopifnot(inherits(dem_rast, "SpatRaster"), inherits(mask_r, "SpatRaster"))
  
  # Align DEM to mask grid (CRS/resolution/extent)
  if (!terra::same.crs(dem_rast, mask_r)) dem_rast <- terra::project(dem_rast, mask_r)
  if (!terra::compareGeom(dem_rast, mask_r, stopOnError = FALSE)) {
    dem_rast <- terra::resample(dem_rast, mask_r, method = "bilinear")
  }
  
  dem_rast <- terra::mask(dem_rast, mask_r)
  
  # Slope magnitude (rise/run). terrain() gives slope angle; tan(angle) => gradient
  slope_rad  <- terra::terrain(dem_rast, v = "slope", unit = "radians", neighbors = 8)
  slope_grad <- base::tan(slope_rad)
  s <- abs(slope_grad)
  
  # Tobler speed (km/h)
  v_kmh <- v0_kmh * exp(-a * abs(s + b))
  v_kmh <- terra::clamp(v_kmh, lower = min_speed_kmh, values = TRUE)
  
  # Convert to resistance (seconds per metre)
  v_ms <- v_kmh * (1000 / 3600)
  resistance <- 1 / v_ms
  
  resistance <- terra::mask(resistance, mask_r)
  resistance
}