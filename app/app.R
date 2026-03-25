library(shiny)
library(sf)
library(terra)
library(weightedVoronoi)

crs_use <- 32636

make_boundary <- function(type = c("square", "concave", "corridor", "rounded_irregular", "split"),
                          crs = crs_use) {
  type <- match.arg(type)
  
  if (type == "rounded_irregular") {
    theta <- seq(0, 2 * pi, length.out = 120)
    r <- 450 +
      70 * sin(3 * theta) +
      45 * cos(5 * theta) +
      30 * sin(7 * theta)
    
    x <- 600 + r * cos(theta)
    y <- 450 + 0.75 * r * sin(theta)
    
    coords <- cbind(x, y)
    coords <- rbind(coords, coords[1, , drop = FALSE])
    
    return(
      st_sf(
        geometry = st_sfc(st_polygon(list(coords))),
        crs = crs
      )
    )
  }
  
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
  
  st_sf(
    geometry = st_sfc(st_polygon(list(coords))),
    crs = crs
  )
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
  plot(
    st_geometry(boundary_sf),
    col = "grey95",
    border = "black",
    lwd = 2,
    main = title
  )
  plot(out$polygons["generator_id"], add = TRUE, border = NA)
  plot(st_geometry(boundary_sf), add = TRUE, border = "black", lwd = 2)
  plot(st_geometry(points_sf), add = TRUE, pch = 21, bg = "red", cex = 1.2)
  text(st_coordinates(points_sf), labels = points_sf$population, pos = 3, cex = 0.85)
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
        "Rounded irregular" = "rounded_irregular",
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
        tabPanel(
          "Map",
          plotOutput("map_plot", height = 550),
          br(),
          downloadButton("download_png", "Download PNG")
        ),
        tabPanel("Summary", tableOutput("summary_tbl")),
        tabPanel("Diagnostics", tableOutput("diagnostics_tbl")),
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
            weight_transform = log10,
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
    req(input$method)
    
    switch(
      input$method,
      
      euclidean = paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: euclidean",
        "Best for: straight-line influence in unconstrained space",
        sep = "\n"
      ),
      
      geodesic = paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic",
        "Engine: classic by default; multisource for additive isotropic repeated runs",
        "Best for: domains where shape constrains access or interaction",
        sep = "\n"
      ),
      
      resistance = paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with resistance_rast",
        "Best for: land cover, friction, or infrastructure effects",
        sep = "\n"
      ),
      
      dem = paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with dem_rast and use_tobler = TRUE",
        "Best for: isotropic slope-dependent movement cost",
        sep = "\n"
      ),
      
      anisotropy = paste(
        "Recommended function: weighted_voronoi_domain()",
        "Distance: geodesic with anisotropy = 'terrain'",
        "Best for: directional uphill/downhill asymmetry",
        sep = "\n"
      ),
      
      "No recommendation available."
    )
  })
  
  output$code_snippet <- renderText({
    req(input$method)
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
        "  weight_model = \"multiplicative\"\n",
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
        "  distance = \"geodesic\"\n",
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
        "  resistance_rast = resistance_rast\n",
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
        "  use_tobler = TRUE\n",
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
        "  downhill_factor = 1.2\n",
        ")"
      ),
      "No code snippet available."
    )
  })
  
  output$diagnostics_tbl <- renderTable({
    s <- scenario()
    
    d <- s$out$diagnostics
    if (is.null(d)) return(NULL)
    
    data.frame(
      metric = names(d),
      value = vapply(d, function(x) {
        if (length(x) == 1) {
          if (is.logical(x)) return(as.character(x))
          if (is.numeric(x)) return(format(signif(x, 4), trim = TRUE))
          return(as.character(x))
        }
        paste0("<", class(x)[1], ">")
      }, character(1)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  
  output$download_png <- downloadHandler(
    filename = function() {
      paste0("weightedVoronoi_", input$method, "_", input$boundary_type, ".png")
    },
    content = function(file) {
      s <- scenario()
      png(file, width = 1400, height = 1000, res = 160)
      plot_tessellation(
        boundary_sf = s$boundary_sf,
        points_sf = s$points_sf,
        out = s$out,
        title = sprintf("%s (%.2fs)", s$method_label, s$elapsed)
      )
      dev.off()
    }
  )
}

shinyApp(ui, server)
