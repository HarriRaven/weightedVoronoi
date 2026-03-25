library(shiny)
library(sf)
library(terra)
library(weightedVoronoi)

crs_use <- 32636

make_boundary <- function(type = c("square", "concave", "corridor", "split"), crs = crs_use) {
  type <- match.arg(type)
  
  coords <- switch(
    type,
    square = rbind(
      c(0, 0), c(1200, 0), c(1200, 900), c(0, 900), c(0, 0)
    ),
    concave = rbind(
      c(0, 0), c(1200, 0), c(1200, 900), c(800, 900), c(800, 500),
      c(400, 500), c(400, 900), c(0, 900), c(0, 0)
    ),
    corridor = rbind(
      c(0, 0), c(500, 0), c(500, 300), c(700, 300), c(700, 0),
      c(1200, 0), c(1200, 900), c(700, 900), c(700, 600),
      c(500, 600), c(500, 900), c(0, 900), c(0, 0)
    ),
    split = rbind(
      c(0, 0), c(500, 0), c(500, 350), c(700, 350), c(700, 0),
      c(1200, 0), c(1200, 900), c(700, 900), c(700, 550),
      c(500, 550), c(500, 900), c(0, 900), c(0, 0)
    )
  )
  
  st_sf(geometry = st_sfc(st_polygon(list(coords))), crs = crs)
}

make_points <- function(boundary_sf, n_points, seed, weight_spread = c("low", "medium", "high")) {
  weight_spread <- match.arg(weight_spread)
  set.seed(seed)
  
  pts <- st_sample(boundary_sf, size = n_points, type = "random")
  
  sdlog <- switch(weight_spread, low = 0.2, medium = 0.5, high = 0.9)
  pop <- round(exp(rnorm(length(pts), log(100), sdlog)))
  pop[pop < 1] <- 1
  
  st_sf(
    id = seq_along(pts),
    population = pop,
    geometry = pts,
    crs = st_crs(boundary_sf)
  )
}

make_demo_resistance <- function(template_rast) {
  r <- template_rast
  terra::values(r) <- 1
  xy <- terra::xyFromCell(r, 1:terra::ncell(r))
  band <- xy[, 1] > 480 & xy[, 1] < 720
  vals <- terra::values(r)
  vals[band] <- 8
  terra::values(r) <- vals
  r
}

make_demo_dem <- function(boundary_sf, res) {
  bnd_v <- terra::vect(boundary_sf)
  dem <- terra::rast(ext = terra::ext(bnd_v), resolution = res, crs = terra::crs(bnd_v))
  xy <- terra::crds(dem, df = TRUE)
  y0 <- (min(xy$y) + max(xy$y)) / 2
  sigma <- (max(xy$y) - min(xy$y)) / 7
  height <- 1000
  terra::values(dem) <- height * exp(-((xy$y - y0)^2) / (2 * sigma^2)) + xy$x * 0.8
  dem
}

plot_tessellation <- function(boundary_sf, points_sf, out, title = "Tessellation") {
  plot(st_geometry(boundary_sf), border = "black", lwd = 2, main = title)
  if (!is.null(out$polygons) && nrow(out$polygons) > 0) {
    plot(out$polygons["generator_id"], add = TRUE)
  }
  plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red", cex = 1.1)
  text(st_coordinates(points_sf), labels = points_sf$population, pos = 3, cex = 0.8)
}

ui <- fluidPage(
  titlePanel("weightedVoronoi demo app (v1)"),
  sidebarLayout(
    sidebarPanel(
      h4("Scenario setup"),
      selectInput("boundary_type", "Boundary", c(
        "Square" = "square",
        "Concave" = "concave",
        "Corridor" = "corridor",
        "Split bottleneck" = "split"
      )),
      sliderInput("n_points", "Number of points", min = 1, max = 10, value = 5, step = 1),
      numericInput("seed", "Random seed", value = 42, min = 1, step = 1),
      selectInput("weight_spread", "Weight spread", c("low", "medium", "high"), selected = "medium"),
      sliderInput("res", "Resolution", min = 10, max = 60, value = 20, step = 5),
      checkboxInput("advanced", "Show advanced methods", value = FALSE),
      conditionalPanel(
        condition = "input.advanced == true",
        radioButtons("method", "Method", choices = c(
          "Euclidean" = "euclidean",
          "Geodesic" = "geodesic",
          "Geodesic + resistance" = "resistance",
          "Geodesic + DEM" = "dem",
          "Terrain anisotropy" = "anisotropy"
        ), selected = "geodesic")
      ),
      conditionalPanel(
        condition = "input.advanced == false",
        radioButtons("method", "Method", choices = c(
          "Euclidean" = "euclidean",
          "Geodesic" = "geodesic"
        ), selected = "geodesic")
      ),
      actionButton("rerun", "Generate example"),
      hr(),
      h4("Recommended workflow"),
      verbatimTextOutput("recommendation")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Map", plotOutput("map_plot", height = 550)),
        tabPanel("Summary", tableOutput("summary_tbl")),
        tabPanel("Code", verbatimTextOutput("code_snippet"))
      )
    )
  )
)

server <- function(input, output, session) {
  scenario <- eventReactive(input$rerun, {
    boundary_sf <- make_boundary(input$boundary_type)
    points_sf <- make_points(boundary_sf, input$n_points, input$seed, input$weight_spread)
    
    base_args <- list(
      points_sf = points_sf,
      weight_col = "population",
      boundary_sf = boundary_sf,
      res = input$res,
      verbose = FALSE
    )
    
    method_label <- switch(
      input$method,
      euclidean = "Euclidean",
      geodesic = "Geodesic",
      resistance = "Geodesic + resistance",
      dem = "Geodesic + DEM",
      anisotropy = "Terrain anisotropy"
    )
    
    runtime <- system.time({
      out <- switch(
        input$method,
        euclidean = do.call(weighted_voronoi_domain, c(base_args, list(
          distance = "euclidean",
          weight_transform = log10,
          weight_model = "multiplicative"
        ))),
        geodesic = do.call(weighted_voronoi_domain, c(base_args, list(
          distance = "geodesic",
          weight_transform = log10,
          geodesic_engine = "classic",
          close_mask = FALSE
        ))),
        resistance = {
          tmp <- do.call(weighted_voronoi_domain, c(base_args, list(
            distance = "euclidean",
            weight_transform = log10,
            weight_model = "multiplicative"
          )))
          rr <- make_demo_resistance(tmp$allocation)
          do.call(weighted_voronoi_domain, c(base_args, list(
            distance = "geodesic",
            weight_transform = log10,
            resistance_rast = rr,
            geodesic_engine = "classic",
            close_mask = FALSE
          )))
        },
        dem = {
          dem_rast <- make_demo_dem(boundary_sf, input$res)
          do.call(weighted_voronoi_domain, c(base_args, list(
            distance = "geodesic",
            weight_transform = log10,
            dem_rast = dem_rast,
            use_tobler = TRUE,
            geodesic_engine = "classic",
            close_mask = FALSE
          )))
        },
        anisotropy = {
          dem_rast <- make_demo_dem(boundary_sf, input$res)
          do.call(weighted_voronoi_domain, c(base_args, list(
            distance = "geodesic",
            dem_rast = dem_rast,
            use_tobler = TRUE,
            anisotropy = "terrain",
            uphill_factor = 3,
            downhill_factor = 1.2,
            geodesic_engine = "classic",
            close_mask = FALSE
          )))
        }
      )
    })
    
    list(
      boundary_sf = boundary_sf,
      points_sf = points_sf,
      out = out,
      method_label = method_label,
      elapsed = unname(runtime["elapsed"])
    )
  }, ignoreNULL = FALSE)
  
  output$map_plot <- renderPlot({
    s <- scenario()
    plot_tessellation(
      boundary_sf = s$boundary_sf,
      points_sf = s$points_sf,
      out = s$out,
      title = sprintf("%s (%.2fs)", s$method_label, s$elapsed)
    )
  })
  
  output$summary_tbl <- renderTable({
    s <- scenario()
    s$out$summary
  })
  
  output$recommendation <- renderText({
    method <- input$method
    if (method == "euclidean") {
      paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: euclidean",
        "Best for: straight-line influence in unconstrained space",
        sep = "\n"
      )
    } else if (method == "geodesic") {
      paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic",
        "Engine: classic for general use; multisource for additive isotropic repeated runs",,
        "Best for: domains where shape constrains access or interaction",
        sep = "\n"
      )
    } else if (method == "resistance") {
      paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with resistance_rast",
        "Best for: land cover, friction, or infrastructure effects",
        sep = "\n"
      )
    } else if (method == "dem") {
      paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with dem_rast and use_tobler = TRUE",
        "Best for: isotropic slope-dependent movement cost",
        sep = "\n"
      )
    } else {
      paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with anisotropy = 'terrain'",
        "Best for: directional uphill/downhill asymmetry",
        sep = "\n"
      )
    }
  })
  
  output$code_snippet <- renderText({
    s <- scenario()
    pts_txt <- paste0("# Example generated with seed = ", input$seed, "\n")
    
    switch(
      input$method,
      euclidean = paste0(
        pts_txt,
        "weighted_voronoi_domain(\n",
        "  points_sf = points_sf,\n",
        "  weight_col = \"population\",\n",
        "  boundary_sf = boundary_sf,\n",
        "  res = ", input$res, ",\n",
        "  weight_transform = log10,\n",
        "  distance = \"euclidean\",\n",
        "  weight_model = \"multiplicative\",\n",
        "  clip_to_boundary = TRUE\n",
        ")"
      ),
      geodesic = paste0(
        pts_txt,
        "weighted_voronoi_domain(\n",
        "  points_sf = points_sf,\n",
        "  weight_col = \"population\",\n",
        "  boundary_sf = boundary_sf,\n",
        "  res = ", input$res, ",\n",
        "  weight_transform = log10,\n",
        "  distance = \"geodesic\",\n",
        "  close_mask = FALSE\n",
        ")"
      ),
      resistance = paste0(
        pts_txt,
        "weighted_voronoi_domain(\n",
        "  points_sf = points_sf,\n",
        "  weight_col = \"population\",\n",
        "  boundary_sf = boundary_sf,\n",
        "  res = ", input$res, ",\n",
        "  weight_transform = log10,\n",
        "  distance = \"geodesic\",\n",
        "  resistance_rast = resistance_rast,\n",
        "  close_mask = FALSE\n",
        ")"
      ),
      dem = paste0(
        pts_txt,
        "weighted_voronoi_domain(\n",
        "  points_sf = points_sf,\n",
        "  weight_col = \"population\",\n",
        "  boundary_sf = boundary_sf,\n",
        "  res = ", input$res, ",\n",
        "  weight_transform = log10,\n",
        "  distance = \"geodesic\",\n",
        "  dem_rast = dem_rast,\n",
        "  use_tobler = TRUE,\n",
        "  close_mask = FALSE\n",
        ")"
      ),
      anisotropy = paste0(
        pts_txt,
        "weighted_voronoi_domain(\n",
        "  points_sf = points_sf,\n",
        "  weight_col = \"population\",\n",
        "  boundary_sf = boundary_sf,\n",
        "  res = ", input$res, ",\n",
        "  distance = \"geodesic\",\n",
        "  dem_rast = dem_rast,\n",
        "  use_tobler = TRUE,\n",
        "  anisotropy = \"terrain\",\n",
        "  uphill_factor = 3,\n",
        "  downhill_factor = 1.2,\n",
        "  close_mask = FALSE\n",
        ")"
      )
    )
  })
}

shinyApp(ui, server)
