---
title: "weightedVoronoi: Weighted Euclidean and geodesic tessellations for ecological spatial allocation within constrained domains"
tags:
  - R
  - spatial analysis
  - landscape ecology
  - Voronoi tessellation
  - geodesic distance
  - least-cost distance
authors:
  - name: Harri Ravenscroft
    affiliation: 1
affiliations:
  - name: Department of Earth Sciences, University of Oxford, Oxford, United Kingdom
    index: 1
date: 2026
bibliography: paper.bib
---

# Summary

Many ecological and social–ecological analyses require dividing space into zones of influence around point locations such as settlements, sampling sites, nests, or service centres. A common mathematical tool for this purpose is a Voronoi tessellation, which partitions space so that every location is assigned to its nearest generator point. However, standard Voronoi tessellations assume that all generators have equal importance and that distance is measured in straight lines across unconstrained space. These assumptions are often unrealistic in ecological systems where sites differ in size or influence and where coastlines, mountains, lakes, or fragmented habitats constrain movement and access.

`weightedVoronoi` is an R package that generates multiplicatively weighted tessellations within arbitrary polygonal domains under either Euclidean or domain-constrained geodesic distance metrics. Geodesic distances may optionally incorporate resistance surfaces, enabling terrain-aware spatial allocation. The package produces complete, connected tessellations that conform to irregular boundaries and returns diagnostic outputs to support reproducible ecological spatial analysis.

# Statement of Need

Ecological spatial allocation problems commonly involve assigning land cover, environmental variables, or management responsibilities to settlements, communities, or sampling sites. In such contexts, allocation must often account for:

1. heterogeneous generator importance (e.g., population size, sampling effort, service capacity),
2. irregular or concave spatial domains (e.g., protected areas, coastlines, fragmented habitats),
3. non-Euclidean accessibility shaped by barriers or terrain.

Straight-line distance may misrepresent ecological interaction or accessibility in mountainous regions, island systems, or landscapes fragmented by unsuitable habitat. While researchers can compute least-cost or constrained distances, translating these into complete spatial partitions typically requires bespoke workflows.

`weightedVoronoi` was developed to provide a unified and reproducible solution for weighted spatial allocation under Euclidean and domain-constrained geodesic distance assumptions. The target audience includes spatial ecologists, landscape ecologists, conservation planners, and researchers linking spatial covariates to socially or biologically defined influence zones.

# State of the Field

Unweighted Voronoi tessellations are available in R through packages such as `deldir`, and planar tessellation utilities are included in spatial point-pattern frameworks such as `spatstat`. These tools generally assume Euclidean distance and do not integrate multiplicative generator weighting within arbitrary polygonal domains. 

Grid-based shortest-path and least-cost distances can be computed using `gdistance`, but this functionality is not packaged as a complete tessellation framework and does not automatically produce spatial partitions or enforce connectivity constraints. In practice, researchers often combine multiple tools or rely on proprietary geographic information system software to implement weighted or constrained tessellations.

`weightedVoronoi` differs by integrating multiplicative weighting, arbitrary polygonal domain constraints, optional resistance-modified geodesic distance, automated enforcement of single connected regions per generator, and structured diagnostics within a single open-source R workflow. This integration provides functionality not simultaneously available in existing R spatial packages.

# Software Design

The package is implemented in R (version 1.0.0) and builds on open-source spatial infrastructure including `sf` for vector data handling [@pebesma2018sf], `terra` for raster operations [@hijmans2025terra], and optionally `gdistance` for graph-based shortest-path computation on grids [@vanetten2017gdistance].

A raster-based architecture was chosen to allow tessellations within complex, concave, or fragmented domains and to enable geodesic and resistance-aware distance calculations. While analytic Voronoi constructions are efficient in unconstrained planar settings, they are less flexible for arbitrary polygonal masks and non-Euclidean distances. The raster approach provides geometric generality at the cost of resolution-dependent computation time.

The core function, `weighted_voronoi_domain()`, accepts generator points with weights, a polygonal boundary, a raster resolution, a weight transformation function, and a distance mode. Allocation follows a multiplicative rule in which each raster cell is assigned to the generator minimizing weighted distance. Post-processing enforces single connected regions per generator and fills unreachable cells in geodesic mode to guarantee complete domain coverage.

This design prioritizes transparency, reproducibility, and flexibility for ecological applications where boundary geometry and accessibility constraints are central.

# Research Impact Statement

The package provides fully reproducible workflows, example vignettes, and unit tests to support immediate application in ecological research. It is distributed under the MIT license and archived with a versioned DOI to facilitate citation and reuse.

The software has been developed in the context of ecological spatial allocation problems and is designed to support applications including settlement boundary approximation, landscape-level influence modelling, and terrain-aware accessibility analysis. By integrating weighting, boundary constraints, and geodesic distance within a single workflow, the package reduces reliance on ad hoc scripting or proprietary software and promotes reproducible ecological spatial modelling.

# Example

```r
library(sf)
library(weightedVoronoi)

crs_use <- 32636

boundary_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0), c(1200, 0), c(1200, 800), c(0, 800), c(0, 0)
  )))),
  crs = crs_use
)

set.seed(1)
pts <- st_sample(boundary_sf, size = 8)

points_sf <- st_sf(
  population = round(exp(rnorm(length(pts), log(300), 0.8))),
  geometry = pts,
  crs = crs_use
)

out <- weighted_voronoi_domain(
  points_sf = points_sf,
  weight_col = "population",
  boundary_sf = boundary_sf,
  res = 20,
  weight_transform = log10,
  distance = "euclidean",
  verbose = FALSE
)
```

## Availability

Source code: https://github.com/HarriRaven/weightedVoronoi

Archived release (v1.0.0): https://doi.org/10.5281/zenodo.18492446

License: MIT


## AI Usage Disclosure

Limited generative AI tools were used to assist in drafting and editing the manuscript text. 
All technical descriptions, software functionality, and methodological content were reviewed and verified by the author to ensure correctness and fidelity to the implemented code. 
No AI tools were used in the development of the software itself.


## Acknowledgements

H. Ravenscroft acknowledges support from the Department of Earth Sciences, University of Oxford, and the Royal Botanic Gardens, Kew.