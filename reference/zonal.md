# Zonal statistics

Compute zonal statistics, that is summarize values of a SpatRaster for
each "zone" defined by another SpatRaster, or by a SpatVector with
polygon geometry.

If `fun` is a true R `function`, the \<SpatRaster,SpatRaster\> method
may fail when using very large SpatRasters, except for the functions
("mean", "min", "max", "sum", "isNA", and "notNA").

You can also summarize values of a SpatVector for each polygon (zone)
defined by another SpatVector.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
zonal(x, z, fun="mean", ..., w=NULL, wide=TRUE,
    as.raster=FALSE, filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRaster,SpatVector'
zonal(x, z, fun="mean", na.rm=FALSE, w=NULL, weights=FALSE,
    exact=FALSE, touches=FALSE, small=TRUE, as.raster=FALSE,
    as.polygons=FALSE, wide=TRUE, filename="", wopt=list())

# S4 method for class 'SpatVector,SpatVector'
zonal(x, z, fun=mean, ..., weighted=FALSE, as.polygons=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- z:

  SpatRaster with cell-values representing zones or a SpatVector with
  each polygon geometry representing a zone. `z` can have multiple
  layers to define intersecting zones

- fun:

  function to be applied to summarize the values by zone. Either as
  character: "mean", "min", "max", "sum", "isNA", and "notNA" and, for
  relatively small SpatRasters, a proper function

- ...:

  additional arguments passed to fun, such as `na.rm=TRUE`

- w:

  SpatRaster with weights. Should have a single-layer with non-negative
  values

- wide:

  logical. Should the values returned in a wide format? For the
  `SpatRaster, SpatRaster` method this only affects the results when
  `nlyr(z) == 2`. For the `SpatRaster, SpatVector` method this only
  affects the results when `fun=table`

- as.raster:

  logical. If `TRUE`, a SpatRaster is returned with the zonal statistic
  for each zone

- filename:

  character. Output filename (ignored if `as.raster=FALSE`

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

- weights:

  logical. If `TRUE` and `y` has polygons, the approximate fraction of
  each cell that is covered is returned as well, for example to compute
  a weighted mean

- exact:

  logical. If `TRUE` and `y` has polygons, the exact fraction of each
  cell that is covered is returned as well, for example to compute a
  weighted mean

- touches:

  logical. If `TRUE`, values for all cells touched by lines or polygons
  are extracted, not just those on the line render path, or whose center
  point is within the polygon. Not relevant for points; and always
  considered `TRUE` when `weights=TRUE` or `exact=TRUE`

- small:

  logical. If `TRUE`, values for all cells in touched polygons are
  extracted if none of the cells center points is within the polygon;
  even if `touches=FALSE`

- weighted:

  logical. If `TRUE`, a weighted.mean is computed and `fun` is ignored.
  Weights are based on the length of the lines or the area of the
  polygons in `x` that intersect with `z`. This argument is ignored of
  `x` is a SpatVector or points

- as.polygons:

  logical. Should the zonal statistics be combined with the geometry of
  `z`?

- na.rm:

  logical. If `TRUE`, `NA`s are removed

## Value

A `data.frame` with a value for each zone, or a SpatRaster, or
SpatVector of polygons.

## See also

See [`global`](https://rspatial.github.io/terra/reference/global.md) for
"global" statistics (i.e., all of `x` is considered a single zone),
[`app`](https://rspatial.github.io/terra/reference/app.md) for local
statistics, and
[`extract`](https://rspatial.github.io/terra/reference/extract.md) for
an alternative way to summarize values of a SpatRaster with a
SpatVector. With
[`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md)
you can compute statistics for cell blocks defined by a number of rows
and columns.

## Examples

``` r
### SpatRaster, SpatRaster
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
z <- rast(r)
values(z) <- rep(c(1:2, NA, 3:4), each=20)
names(z) <- "zone"
zonal(r, z, "sum", na.rm=TRUE)
#>   zone lyr.1
#> 1    1   210
#> 2    2   610
#> 3    3  1410
#> 4    4  1810

# with weights 
w <- init(r, "col")
zonal(r, z, w=w, "mean", na.rm=TRUE)
#>   zone lyr.1
#> 1    1    12
#> 2    2    32
#> 3    3    72
#> 4    4    92

# multiple layers
r <- rast(system.file("ex/logo.tif", package = "terra")) 
# zonal layer 
z <- rast(r, 1)
names(z) <- "zone"
values(z) <- rep(c(1:2, NA, c(3:4)), each=ncell(r)/5, length.out=ncell(r))

zonal(r, z, "mean", na.rm = TRUE)
#>   zone      red    green     blue
#> 1    1 197.9486 198.0103 193.5556
#> 2    2 173.2219 176.7717 185.2585
#> 3    3 168.2952 172.6232 184.6939
#> 4    4 193.5859 197.0019 206.5717

# raster of zonal values
zr <- zonal(r, z, "mean", na.rm = TRUE, as.raster=TRUE)


### SpatRaster, SpatVector
x <- rast(ncol=2,nrow=2, vals=1:4, xmin=0, xmax=1, ymin=0, ymax=1, crs="+proj=utm +zone=1")
p <- as.polygons(x)
pp <- shift(p, .2)
r <- disagg(x, 4)

zonal(r, p)
#>   lyr.1
#> 1     1
#> 2     2
#> 3     3
#> 4     4
zonal(r, p, sum)
#>   lyr.1
#> 1    16
#> 2    32
#> 3    48
#> 4    64
zonal(x, pp, exact=TRUE)
#>   lyr.1
#> 1   1.4
#> 2   2.0
#> 3   3.4
#> 4   4.0
zonal(c(x, x*10), pp, w=x)
#>   lyr.1 lyr.1
#> 1     1    10
#> 2     2    20
#> 3     3    30
#> 4     4    40


### SpatVector, SpatVector

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)[,c(2,4)]

p <- spatSample(v, 100)
values(p) <- data.frame(b2=1:100, ssep1=100:1)

zonal(p, v, mean)
#>    zone       b2    ssep1
#> 1     1 52.08333 48.91667
#> 2     2 66.28571 34.71429
#> 3     3 40.07692 60.92308
#> 4     4 51.33333 49.66667
#> 5     5 62.16667 38.83333
#> 6     6 42.33333 58.66667
#> 7     7 48.00000 53.00000
#> 8     8 48.88889 52.11111
#> 9     9 46.16667 54.83333
#> 10   10 50.00000 51.00000
#> 11   11 53.30000 47.70000
#> 12   12 52.45455 48.54545
```
