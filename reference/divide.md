# Subdivide a raster or polygons

Divide a `SpatRaster` into `n` parts with approximately the same sum of
weights (cell values).

Divides a `SpatVector` of polygons into `n` compact and approximately
equal area parts. The results are not deterministic so you should use
set.seed to be able to reproduce your results. If you get a warning
about non-convergence, you can increase the number of iterations used
with additional argument iter.max

## Usage

``` r
# S4 method for class 'SpatRaster'
divide(x, n=2, start="ns", as.raster=FALSE, na.rm=TRUE)

# S4 method for class 'SpatVector'
divide(x, n=5, w=NULL, alpha=1, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector of polygons

- n:

  numeric. Can be a single positive integer to indicate the number of
  parts (SpatVector) or the number of splits (SpatRaster).

  If `x` is a SpatRaster, it can also be a vector with values -2, -1, 1,
  or 2. Where 1 means one split and 2 means two splits, and the negative
  sign indicates an East-West (vertical) split as opposed to a
  North-South split.

  If `x` is a SpatVector it can be a list with at least one of these
  elements: `horizontal` and `vertical` that specify the proportions of
  the area that splits should cover. This can either be a single
  fraction such as 1/3, or a sequence of fractions in ascending order
  such as `c(1/4, 1/2, 1)`

- start:

  character. To indicate the initial direction of splitting the raster.
  "ns" for North-South (horizontal) or "ew" for East-West (vertical)

- as.raster:

  logical. If `FALSE` a SpatVector is returned. If `FALSE`, a SpatRaster
  is returned. If `NA` a list with a SpatRaster and a SpatVector is
  returned

- na.rm:

  logical. If `TRUE` cells in `x` that are `NA` are not included in the
  output

- w:

  SpatRaster with, for example, environmental data

- alpha:

  numeric. One or two numbers that act as weights for the x and y
  coordinates

- ...:

  additional arguments such as `iter.max` passed on to
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## See also

[`thresh`](https://rspatial.github.io/terra/reference/thresh.md);
[`makeTiles`](https://rspatial.github.io/terra/reference/makeTiles.md)

## Value

SpatVector or SpatRaster, or a list with both

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- divide(r, 3)
# plot(r); lines(x)


f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
d <- divide(v, 3)
dv <- divide(v, list(h=.5))
```
