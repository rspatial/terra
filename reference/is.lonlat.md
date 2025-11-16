# Check for longitude/latitude crs

Test whether a SpatRaster or SpatVector has a longitude/latitude
coordinate reference system (CRS), or perhaps has one. That is, when the
CRS is unknown (`""`) but the x coordinates are within -181 and 181 and
the y coordinates are within -90.1 and 90.1. For a SpatRaster you can
also test if it has a longitude/latitude CRS and it is "global" (covers
all longitudes).

A warning is given if the CRS is missing or if it is specified as
longitude/latitude but the coordinates do not match that.

## Usage

``` r
# S4 method for class 'SpatRaster'
is.lonlat(x, perhaps=FALSE, warn=TRUE, global=FALSE)

# S4 method for class 'SpatVector'
is.lonlat(x, perhaps=FALSE, warn=TRUE)

# S4 method for class 'character'
is.lonlat(x, perhaps=FALSE, warn=TRUE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- perhaps:

  logical. If `TRUE` and the CRS is unknown, the method returns `TRUE`
  if the coordinates are plausible for longitude/latitude

- warn:

  logical. If `TRUE`, a warning is given if the CRS is unknown but
  assumed to be lon/lat and `perhaps=TRUE`

- global:

  logical. If `TRUE`, the method tests if the raster covers all
  longitudes (from -180 to 180 degrees) such that the extreme columns
  are in fact adjacent

## Value

logical or NA

## Examples

``` r
r <- rast()
is.lonlat(r)
#> [1] TRUE
is.lonlat(r, global=TRUE)
#> [1] TRUE

crs(r) <- ""
is.lonlat(r)
#> Warning: [is.lonlat] unknown crs
#> [1] NA
is.lonlat(r, perhaps=TRUE, warn=FALSE)
#> [1] TRUE

crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
is.lonlat(r)
#> [1] FALSE
```
