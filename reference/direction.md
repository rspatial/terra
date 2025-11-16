# Direction

The direction (azimuth) to or from the nearest cell that is not `NA`.
The direction is expressed in radians, unless you use argument
`degrees=TRUE`.

## Usage

``` r
# S4 method for class 'SpatRaster'
direction(x, from=FALSE, degrees=FALSE, method="cosine", filename="", ...)
```

## Arguments

- x:

  SpatRaster

- from:

  Logical. Default is `FALSE`. If `TRUE`, the direction from (instead of
  to) the nearest cell that is not `NA` is returned

- degrees:

  Logical. If `FALSE` (the default) the unit of direction is radians.

- method:

  character. Should be "geo", or "cosine". With "geo" the most precise
  but slower geodesic method of Karney (2003) is used. The "cosine"
  method is faster but less precise

- filename:

  Character. Output filename (optional)

- ...:

  Additional arguments as for
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`distance`](https://rspatial.github.io/terra/reference/distance.md)

## Examples

``` r
r <- rast(ncol=36,nrow=18, crs="+proj=merc")
values(r) <- NA
r[306] <- 1
b <- direction(r, degrees=TRUE) 
plot(b)


crs(r) <- "+proj=longlat"
b <- direction(r) 
plot(b)

```
