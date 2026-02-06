# extract values along lines

Extract raster values along a line. That is, the returned values are
ordered along the line. That is not the case with
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Usage

``` r
extractAlong(x, y, ID=TRUE, cells=FALSE, xy=FALSE, online=FALSE, bilinear=TRUE)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector with lines geometry

- ID:

  logical. Should an ID column be added? If so, the first column
  returned has the IDs (record numbers) of input SpatVector `y`

- cells:

  logical. If `TRUE` the cell numbers are also returned

- xy:

  logical. If `TRUE` the coordinates of the cells traversed by `y` are
  also returned. See
  [`xyFromCell`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

- online:

  logical. If `TRUE` the returned coordinates are snapped to `y`

- bilinear:

  logical. If `TRUE` the returned raster values computed with bilinear
  interpolation from the nearest four cells. Only relevant if
  `online=TRUE`

## Value

data.frame

## See also

[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
r <- rast(ncols=36, nrows=18, vals=1:(18*36))
cds1 <- rbind(c(-50,0), c(0,60), c(40,5), c(15,-45), c(-10,-25))
cds2 <- rbind(c(80,20), c(140,60), c(160,0), c(140,-55))
lines <- vect(list(cds1, cds2), "lines")

extractAlong(r, lines)
#> Warning: [is.lonlat] unknown crs
#> Error in if (is.lonlat(v, warn = FALSE) && is.lonlat(x, warn = FALSE)) {    crs(v) <- NULL}: missing value where TRUE/FALSE needed
```
