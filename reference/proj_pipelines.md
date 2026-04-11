# Find CRS transformation pipelines

Find candidate coordinate transformation pipelines between two
coordinate reference systems. This can be used to examine which
transformations are available, whether they require grid files, and what
accuracy they provide.

## Usage

``` r
proj_pipelines(from, to, authority="", AOI=NULL, use="NONE", grid_availability="USED", 
  desired_accuracy=-1.0, strict_containment=FALSE, axis_order_authority_compliant=FALSE)
```

## Arguments

- from:

  character, SpatRaster, SpatVector, or any object from which
  [`crs()`](https://rspatial.github.io/terra/reference/crs.md) can
  extract a CRS

- to:

  same types as for `from`

- authority:

  character. Constrain output pipelines to those from this authority
  (e.g. `"EPSG"`). Default `""` means no constraint

- AOI:

  SpatExtent. Area of interest for the transformation in degrees. For an
  area of interest crossing the anti-meridian, west should be greater
  than east, e.g. `ext(150, 210, 0, 30)` usecharacter. One of `"NONE"`,
  `"BOTH"`, `"INTERSECTION"`, `"SMALLEST"`, indicating how areas of
  interest of source and target CRS are used. Ignored when `AOI` is
  specified grid_availabilitycharacter. One of `"USED"` (grid
  availability is used for sorting; operations with missing grids are
  sorted last), `"DISCARD"` (discard operations if a required grid is
  missing), `"IGNORED"` (ignore grid availability; results presented as
  if all grids were available), or `"AVAILABLE"` (results presented as
  if grids known to PROJ were available; typical when networking is
  enabled) desired_accuracynumeric. Only return pipelines with at least
  this accuracy (in metres). Use `-1` (the default) to impose no
  accuracy constraint strict_containmentlogical. If `FALSE` (the
  default), permit partial matching of the area of interest. If `TRUE`,
  the area of interest must be strictly contained
  axis_order_authority_compliantlogical. If `FALSE` (the default),
  always use longitude/easting for the first axis. If `TRUE`, follow the
  axis order given by the CRS definitions

A `data.frame` with "pipelines" and the following columns:

- description:

  human-readable description

- definition:

  PROJ pipeline definition string

- has_inverse:

  logical; whether the pipeline has an inverse

- accuracy:

  accuracy in metres; `-1` indicates unknown (ballpark) accuracy

- grid_count:

  number of grids required by this pipeline

- instantiable:

  logical; whether the pipeline can be used (all required grids are
  available)

Requires PROJ \>= 7.1 (inc" 3). See
[`projNetwork`](https://rspatial.github.io/terra/reference/gdal.md) to
enable online grid downloading.

[`crs`](https://rspatial.github.io/terra/reference/crs.md),
[`project`](https://rspatial.github.io/terra/reference/project.md),
[`same.crs`](https://rspatial.github.io/terra/reference/same.crs.md),
[`projNetwork`](https://rspatial.github.io/terra/reference/gdal.md)

\# Simple case: WGS84 to UTM zone 32Nproj_pipelines("EPSG:4326",
"EPSG:32632")# Many pipelines between NAD83/Conus Albers and WGS84 pp
\<- proj_pipelines("EPSG:5070", "EPSG:4326") nrow(pp)# Only pipelines
that work without downloading grids pp_local \<-
proj_pipelines("EPSG:5070", "EPSG:4326", grid_availability="DISCARD")
pp_local# Using a SpatRaster as input r \<- rast(ncols=10, nrows=10,
crs="EPSG:4326") proj_pipelines(r, "EPSG:3857")

spatial
