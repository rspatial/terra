# Geolocation arrays and GCPs

Detect whether a SpatRaster source has GDAL geolocation arrays (e.g.
satellite swath / curvilinear NetCDF) and/or Ground Control Points
(GCPs), and report the associated spatial reference and array paths.

When a file has no conventional geotransform but does have a
`GEOLOCATION` metadata domain (or GCPs),
[`rast`](https://rspatial.github.io/terra/reference/rast.md) stores that
information on the source and, if the dataset itself has no CRS, applies
the geolocation/GCP CRS automatically so that
[`project`](https://rspatial.github.io/terra/reference/project.md) can
use the arrays.

## Usage

``` r
# S4 method for class 'SpatRaster'
has.geoloc(x)

# S4 method for class 'SpatRaster'
geoloc(x)
```

## Arguments

- x:

  SpatRaster

## Value

`has.geoloc`: logical, one value per raster data source.

`geoloc`: `data.frame` with one row per source and columns `geolocation`
(logical), `gcps` (logical), `srs` (character), `x` / `y` (paths to the
X/Y geolocation datasets when present).

## See also

[`project`](https://rspatial.github.io/terra/reference/project.md)`, `[`rectify`](https://rspatial.github.io/terra/reference/rectify.md)`, `[`is.rotated`](https://rspatial.github.io/terra/reference/is.rotated.md)`, `[`crs`](https://rspatial.github.io/terra/reference/crs.md)

## Examples

``` r
r <- rast(nrows=10, ncols=10)
has.geoloc(r)
#> [1] FALSE
geoloc(r)
#>   geolocation  gcps srs x y
#> 1       FALSE FALSE        
```
