# Find PROJ transformation pipelines

Find candidate coordinate transformation pipelines between two
coordinate reference systems using the PROJ library. This can be used to
examine which transformations are available, whether they require grid
files, and what accuracy they provide.

## Usage

``` r
proj_pipelines(from, to, authority="", AOI=numeric(0), use="NONE", 
  grid_availability="USED", desired_accuracy=-1.0, strict_containment=FALSE,
    axis_order_authority_compliant=FALSE)
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

  numeric of length four. Area of interest for the transformation:
  (west, south, east, north) in degrees. For an area of interest
  crossing the anti-meridian, west will be greater than east

- use:

  character. One of `"NONE"`, `"BOTH"`, `"INTERSECTION"`, `"SMALLEST"`,
  indicating how areas of interest of source and target CRS are used.
  Ignored when `AOI` is specified

- grid_availability:

  character. One of `"USED"` (grid availability is used for sorting;
  operations with missing grids are sorted last), `"DISCARD"` (discard
  operations if a required grid is missing), `"IGNORED"` (ignore grid
  availability; results presented as if all grids were available), or
  `"AVAILABLE"` (results presented as if grids known to PROJ were
  available; typical when networking is enabled)

- desired_accuracy:

  numeric. Only return pipelines with at least this accuracy (in
  metres). Use `-1` (the default) to impose no accuracy constraint

- strict_containment:

  logical. If `FALSE` (the default), permit partial matching of the area
  of interest. If `TRUE`, the area of interest must be strictly
  contained

- axis_order_authority_compliant:

  logical. If `FALSE` (the default), always use longitude/easting for
  the first axis. If `TRUE`, follow the axis order given by the CRS
  definitions

## Value

A `data.frame` with class `"proj_pipelines"` and the following columns:

- id:

  pipeline identifier (may be empty)

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

Attributes `"from"` and `"to"` store the CRS strings used.

## Note

Requires PROJ \>= 7.1 (included with GDAL \>= 3). See
[`projNetwork`](https://rspatial.github.io/terra/reference/gdal.md) to
enable online grid downloading.

## See also

[`crs`](https://rspatial.github.io/terra/reference/crs.md),
[`project`](https://rspatial.github.io/terra/reference/project.md),
[`same.crs`](https://rspatial.github.io/terra/reference/same.crs.md),
[`projNetwork`](https://rspatial.github.io/terra/reference/gdal.md)

## Examples

``` r
# Simple case: WGS84 to UTM zone 32N
proj_pipelines("EPSG:4326", "EPSG:32632")
#>                             description
#> 1 axis order change (2D) + UTM zone 32N
#>                                                                                            definition
#> 1 +proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=utm +zone=32 +ellps=WGS84
#>   has_inverse accuracy grid_count instantiable
#> 1        TRUE        0          0         TRUE

# Many pipelines between NAD83/Conus Albers and WGS84
pp <- proj_pipelines("EPSG:5070", "EPSG:4326")
nrow(pp)
#> [1] 53

# Only pipelines that work without downloading grids
pp_local <- proj_pipelines("EPSG:5070", "EPSG:4326", grid_availability="DISCARD")
pp_local
#>                                                                                          description
#> 1                             Inverse of Conus Albers + NAD83 to WGS 84 (1) + axis order change (2D)
#> 2                             Inverse of Conus Albers + NAD83 to WGS 84 (3) + axis order change (2D)
#> 3                             Inverse of Conus Albers + NAD83 to WGS 84 (2) + axis order change (2D)
#> 4 Inverse of Conus Albers + Ballpark geographic offset from NAD83 to WGS 84 + axis order change (2D)
#>                                                                                                                                                                                                                                                                                                   definition
#> 1                                                                                                                                                 +proj=pipeline +step +inv +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +step +proj=unitconvert +xy_in=rad +xy_out=deg
#> 2 +proj=pipeline +step +inv +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +step +proj=push +v_3 +step +proj=cart +ellps=GRS80 +step +proj=helmert +x=1 +y=1 +z=-1 +step +inv +proj=cart +ellps=WGS84 +step +proj=pop +v_3 +step +proj=unitconvert +xy_in=rad +xy_out=deg
#> 3 +proj=pipeline +step +inv +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +step +proj=push +v_3 +step +proj=cart +ellps=GRS80 +step +proj=helmert +x=-2 +y=0 +z=4 +step +inv +proj=cart +ellps=WGS84 +step +proj=pop +v_3 +step +proj=unitconvert +xy_in=rad +xy_out=deg
#> 4                                                                                                                                                 +proj=pipeline +step +inv +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +step +proj=unitconvert +xy_in=rad +xy_out=deg
#>   has_inverse accuracy grid_count instantiable
#> 1        TRUE        4          0         TRUE
#> 2        TRUE        4          0         TRUE
#> 3        TRUE        8          0         TRUE
#> 4        TRUE       -1          0         TRUE

# Using a SpatRaster as input
r <- rast(ncols=10, nrows=10, crs="EPSG:4326")
proj_pipelines(r, "EPSG:3857")
#>                                                      description
#> 1 axis order change (2D) + Popular Visualisation Pseudo-Mercator
#>                                                                                                                       definition
#> 1 +proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=webmerc +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
#>   has_inverse accuracy grid_count instantiable
#> 1        TRUE        0          0         TRUE
```
