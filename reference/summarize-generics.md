# Summarize

Compute summary statistics for cells, either across layers or between
layers (parallel summary).

The following summary methods are available for SpatRaster:
`any, anyNA, all, allNA, max, min, mean, median, prod, range, stdev, sum, which.min, which.max`.
See [`modal`](https://rspatial.github.io/terra/reference/modal.md) to
compute the mode and
[`app`](https://rspatial.github.io/terra/reference/app.md) to compute
summary statistics that are not included here.

Because generic functions are used, the method applied is chosen based
on the first argument: "`x`". This means that if `r` is a SpatRaster,
`mean(r, 5)` will work, but `mean(5, r)` will not work.

The `mean` method has an argument "trim" that is ignored.

If `pop=TRUE` `stdev` computes the population standard deviation,
computed as:

`f <- function(x) sqrt(sum((x-mean(x))^2) / length(x))`

This is different than the sample standard deviation returned by `sd`
(which uses `n-1` as denominator).

## Usage

``` r
# S4 method for class 'SpatRaster'
min(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
max(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
range(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
prod(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
sum(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
any(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
all(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
range(x, ..., na.rm=FALSE)

# S4 method for class 'SpatRaster'
which.min(x)

# S4 method for class 'SpatRaster'
which.max(x)

# S4 method for class 'SpatRaster'
stdev(x, ..., pop=TRUE, na.rm=FALSE)

# S4 method for class 'SpatRaster'
mean(x, ..., trim=NA, na.rm=FALSE)

# S4 method for class 'SpatRaster'
median(x, na.rm=FALSE, ...)

# S4 method for class 'SpatRaster'
anyNA(x)

# S4 method for class 'SpatRaster'
countNA(x, n=0)

# S4 method for class 'SpatRaster'
noNA(x, falseNA=FALSE)

# S4 method for class 'SpatRaster'
allNA(x, falseNA=FALSE)
```

## Arguments

- x:

  SpatRaster

- ...:

  additional SpatRasters or numeric values; and arguments `par` for
  parallel summarization (see Details), and `filename`, `overwrite` and
  `wopt` as for
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

- na.rm:

  logical. If `TRUE`, `NA` values are ignored. If `FALSE`, `NA` is
  returned if `x` has any `NA` values

- trim:

  ignored

- pop:

  logical. If `TRUE`, the population standard deviation is computed.
  Otherwise the sample standard deviation is computed

- falseNA:

  logical. If `TRUE`, cells that would otherwise be `FALSE` are set to
  `NA`

- n:

  integer. If `n > 0`, cell values are `TRUE` if at least `n` of its
  layers are `NA`

## Value

SpatRaster

## Details

Additional argument `par` can be used for "parallel" summarizing a
SpatRaster and a numeric or logical value. If a SpatRaster `x` has three
layers, `max(x, 5)` will return a single layer (the number five is
treated as a layer in which all cells have value five). In contrast
`max(x, 5, par=TRUE)` returns three layers (the number five is treated
as another SpatRaster with a single layer with all cells having the
value five.

## See also

[`app`](https://rspatial.github.io/terra/reference/app.md),
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md),
[`modal`](https://rspatial.github.io/terra/reference/modal.md),
[`which.lyr`](https://rspatial.github.io/terra/reference/which.md)

## Examples

``` r
set.seed(0)
r <- rast(nrows=10, ncols=10, nlyrs=3)
values(r) <- runif(ncell(r) * nlyr(r))

x <- mean(r)
# note how this returns one layer
x <- sum(c(r, r[[2]]), 5)

# and this returns three layers
y <- sum(r, r[[2]], 5)

max(r)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :       max 
#> min value   : 0.1808664 
#> max value   : 0.9926841 

## when adding a number, do you want 1 layer or all layers?
# 1 layer
max(r, 0.5)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :       max 
#> min value   : 0.5000000 
#> max value   : 0.9926841 

# all layers
max(r, 0.5, par=TRUE)
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       :     lyr.1,     lyr.2,     lyr.3 
#> min values  : 0.5000000, 0.5000000, 0.5000000 
#> max values  : 0.9919061, 0.9926841, 0.9815635 

y <- stdev(r)
# not the same as 
yy <- app(r, sd)

z <- stdev(r, r*2)

x <- mean(r, filename=paste0(tempfile(), ".tif"))


v <- values(r)
set.seed(3)
v[sample(length(v), 50)] <- NA
values(r) <- v
is.na(r)
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.2, lyr.3 
#> min values  : FALSE, FALSE, FALSE 
#> max values  :  TRUE,  TRUE,  TRUE 
anyNA(r)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :  lyr1 
#> min value   : FALSE 
#> max value   :  TRUE 
allNA(r)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :  lyr1 
#> min value   : FALSE 
#> max value   :  TRUE 
countNA(r)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : lyr1 
#> min value   :    0 
#> max value   :    3 
countNA(r, 2)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :  lyr1 
#> min value   : FALSE 
#> max value   :  TRUE 
```
