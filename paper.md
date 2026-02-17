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

Spatial allocation around point locations is a recurring requirement in ecology and social–ecological research, including delineating settlement boundaries, defining service or management catchments, linking land cover to communities or sampling sites, and modelling influence zones across heterogeneous landscapes. Standard Voronoi (Thiessen) tessellations assume equal generator influence and unconstrained Euclidean space, assumptions that are frequently violated in ecological systems where sites differ in importance and effective distance is shaped by complex boundaries, barriers, and terrain.

`weightedVoronoi` is an R package that generates **multiplicatively weighted tessellations** within **arbitrary polygonal domains** under either **Euclidean** or **domain-constrained geodesic** distance metrics, with optional **resistance-modified geodesic distances** (e.g., terrain-aware allocation). Tessellations are computed on a rasterised representation of the domain to ensure complete coverage of irregular and concave geometries. Post-processing enforces a single connected region per generator and removes enclave artefacts. The package returns polygon and raster tessellations alongside generator-level summaries and diagnostics to support reproducible ecological spatial allocation workflows.

# Statement of Need

Ecological spatial analyses often require assigning each location within a bounded region to the “nearest” or most influential site while accounting for heterogeneous importance (e.g., population size, sampling effort, service capacity) and constrained accessibility. In coastal systems, fragmented habitats, mountainous terrain, or protected areas, Euclidean distance may poorly represent interaction or access. 

Existing R tools address parts of this problem but do not provide an integrated workflow. Unweighted Voronoi tessellations are available via packages such as `deldir`, while spatial point-pattern frameworks such as `spatstat` provide related planar tessellation utilities. Least-cost and grid-based distance calculations can be performed using `gdistance`. However, these approaches typically:

- assume Euclidean distance or unconstrained planar space,
- do not integrate multiplicative generator weighting within arbitrary polygonal domains,
- require bespoke workflows to derive spatial partitions from cost surfaces,
- do not automatically enforce single connected regions per generator,
- lack structured diagnostics describing allocation completeness and connectivity.

`weightedVoronoi` fills this gap by providing a unified, open-source implementation for producing **weighted Euclidean** and **boundary-constrained geodesic** tessellations within arbitrary domains, optionally incorporating resistance surfaces, and returning reproducible diagnostics suitable for ecological sensitivity analysis.

# Software Description

The package (version 1.0.0) is implemented in R and depends on open-source spatial infrastructure including `sf` for vector data handling [@pebesma2018sf], `terra` for raster operations [@hijmans2025terra], and optionally `gdistance` for shortest-path distance computation on grids [@vanetten2017gdistance].

The core function is:

weighted_voronoi_domain()


## Inputs

Users supply:

- an `sf` point object containing generator locations and a numeric weight attribute,
- an `sf` polygon (`boundary_sf`) defining the tessellation domain,
- a raster resolution controlling discretisation,
- a weight transformation function (e.g., identity, logarithmic, power),
- a distance mode (`"euclidean"` or `"geodesic"`).

Inputs must use a projected coordinate reference system with metric units.

## Weighted allocation rule

Each raster cell location \(x\) within the domain \(D\) is assigned to the generator \(p_i\) that minimises:

\[
\arg\min_i \frac{d(x, p_i)}{g(w_i)}
\]

where \(d(\cdot)\) is Euclidean or geodesic distance, \(w_i\) is the generator-specific weight, and \(g(\cdot)\) is a user-defined transformation. This multiplicative formulation supports ecologically motivated assumptions such as diminishing marginal influence via logarithmic or sub-linear scaling.

## Distance modes

All tessellations are computed on a rasterised representation of the domain to ensure exact conformance to irregular boundaries.

- **Euclidean mode**: Straight-line distances are computed efficiently on the raster grid using `terra`.
- **Geodesic mode**: Shortest-path distances are computed within the rasterised domain using graph-based methods via `gdistance`. Cells unreachable from a given generator are assigned infinite cost.

Optionally, geodesic distances may incorporate a resistance surface (e.g., elevation- or slope-derived movement cost) to produce terrain-aware tessellations.

## Post-processing and outputs

To improve interpretability and robustness:

- each generator is guaranteed ownership of the raster cell containing its location,
- disconnected “island” regions are removed to enforce a single connected component per generator,
- cells unreachable from all generators (geodesic mode) are reassigned locally to preserve complete domain coverage.

The function returns:

1. an `sf` polygon tessellation (one polygon per generator),
2. an allocation raster (`terra::SpatRaster`),
3. a generator-level summary table (allocated area, area share, cell count, raw and transformed weights),
4. a structured diagnostics object (domain closure, raster resolution, number of generators retained, unreachable fractions in geodesic mode, and parameter settings).

Given identical inputs and resolution, tessellations are deterministic and reproducible.

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

Availability

Source code: https://github.com/HarriRaven/weightedVoronoi

Archived release (v1.0.0): https://doi.org/10.5281/zenodo.18492446

License: MIT

Acknowledgements

H. Ravenscroft acknowledges support from the Department of Earth Sciences, University of Oxford, and the Royal Botanic Gardens, Kew.