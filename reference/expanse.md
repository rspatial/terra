# Get the expanse (area) of individual polygons or for all (summed) raster cells

Compute the area covered by polygons or for all raster cells that are
not `NA`.

This method computes areas for longitude/latitude rasters, as the size
of the cells is constant in degrees, but not in square meters. But it
can also be important if the coordinate reference system is planar, but
not equal-area.

For vector data, the best way to compute area is to use the
longitude/latitude CRS. This is contrary to (erroneous, but popular)
belief that you should use a planar coordinate reference system. Where
applicable, the transformation to lon/lat is done automatically, if
`transform=TRUE`.

Note that it is important that polygon geometries are valid. If they are
not valid, the computed area may be wrong. You can check for validity
with
[`is.valid`](https://rspatial.github.io/terra/reference/is.valid.md) and
fix some problems with
[`makeValid`](https://rspatial.github.io/terra/reference/is.valid.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
expanse(x, unit="m", transform=TRUE, byValue=FALSE,
      zones=NULL, wide=FALSE, usenames=FALSE)

# S4 method for class 'SpatVector'
expanse(x, unit="m", transform=TRUE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- unit:

  character. Output unit of area. One of "m", "km", or "ha"

- transform:

  logical. If `TRUE`, planar CRS are transformed to lon/lat for accuracy

- byValue:

  logical. If `TRUE`, the area for each unique cell value is returned

- zones:

  NULL or SpatRaster with the same geometry identifying zones in `x`

- wide:

  logical. Should the results be in "wide" rather than "long" format?

- usenames:

  logical. If `TRUE` layers are identified by their names instead of
  their numbers

## Value

**SpatRaster**: `data.frame` with at least two columns ("layer" and
"area") and possibly also "value" (if `byValue` is `TRUE`), and "zone"
(if `zones` is `TRUE`). If `x` has no values, the total area of all
cells is returned. Otherwise, the area of all cells that are not `NA` is
returned.

**SpatVector**: numeric (one value for each (multi-) polygon geometry.

## See also

[`cellSize`](https://rspatial.github.io/terra/reference/cellSize.md) for
a the size of individual cells of a raster, that can be summed with
[`global`](https://rspatial.github.io/terra/reference/global.md) or with
[`zonal`](https://rspatial.github.io/terra/reference/zonal.md) to get
the area for different zones;
[`surfArea`](https://rspatial.github.io/terra/reference/surfArea.md) for
a raster with elevation values, taking into account the sloping nature
of the surface.

## Examples

``` r
### SpatRaster 
r <- rast(nrows=18, ncols=36)
v <- 1:ncell(r)
v[200:400] <- NA
values(r) <- v

# summed area in km2
expanse(r, unit="km")
#>   layer      area
#> 1     1 273986501

# all cells 
expanse(rast(r), unit="km")
#>   layer      area
#> 1     1 510065622

r <- rast(ncols=90, nrows=45, ymin=-80, ymax=80)
m <- project(r, "+proj=merc")

expanse(m, unit="km")
#>   layer      area
#> 1     1 498751903
expanse(m, unit="km", transform=FALSE)
#>   layer       area
#> 1     1 1241591858

m2 <- c(m, m)
values(m2) <- cbind(c(1,2,NA,NA), c(11:14))
#> Warning: [setValues] values were recycled
expanse(m2, unit="km", byValue=TRUE, wide=TRUE)
#>   layer        1        2       11       12       13       14
#> 1     1 10182145 10182145        0        0        0        0
#> 3     2        0        0 10182145 10182145 10182145 10182145


v <- vect(system.file("ex/lux.shp", package="terra"))
r <- rast(system.file("ex/elev.tif", package="terra"))
r <- round((r-50)/100)
levels(r) <- data.frame(id=1:5, name=c("forest", "water", "urban", "crops", "grass"))
expanse(r, byValue=TRUE)
#>   layer  value       area
#> 1     1 forest   50778375
#> 2     1  water  715780779
#> 3     1  urban 1118216460
#> 4     1  crops  617965652
#> 5     1  grass   60868836

g <- rasterize(v, r, "NAME_1")
expanse(r, byValue=TRUE, zones=g, wide=TRUE)
#>   layer         zone   forest     water     urban     crops    grass
#> 1     1     Diekirch  1110530 146734320 332045014 584531527 58106847
#> 2     1 Grevenmacher 46327739 247619382 219364571   1670053        0
#> 5     1   Luxembourg        0 315854211 555664543  25106417        0


### SpatVector
v <- vect(system.file("ex/lux.shp", package="terra"))

a <- expanse(v)
a
#>  [1] 312283206 218674025 259454806  76200409 263174257 188282143 128991500
#>  [8] 210354494 185630770 251322021 237113004 233329960
sum(a)
#> [1] 2564810595
```
