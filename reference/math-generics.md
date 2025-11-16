# General mathematical methods

Standard mathematical methods for computations with SpatRasters.
Computations are local (applied on a cell by cell basis). If multiple
SpatRasters are used, these must have the same extent and resolution.
These have been implemented:

`abs, sign, sqrt, ceiling, floor, trunc, cummax, cummin, cumprod, cumsum, log, log10, log2, log1p, acos, acosh, asin, asinh, atan, atanh, exp, expm1, cos, cosh, sin, sinh, tan, tanh, round, signif`

Instead of directly calling these methods, you can also provide their
name to the `math` method. This is useful if you want to provide an
output filename.

The following methods have been implemented for `SpatExtent`:
`round, floor, ceiling`

`round` has also been implemented for `SpatVector`, to round the
coordinates of the geometries.

## Usage

``` r
# S4 method for class 'SpatRaster'
sqrt(x)

# S4 method for class 'SpatRaster'
log(x, base=exp(1))

# S4 method for class 'SpatRaster'
round(x, digits=0)

# S4 method for class 'SpatRaster'
math(x, fun, digits=0, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatVector'
round(x, digits=4)

# S4 method for class 'SpatRaster'
cumsum(x)
```

## See also

See [`app`](https://rspatial.github.io/terra/reference/app.md) to use
mathematical functions not implemented by the package, and
[`Arith-methods`](https://rspatial.github.io/terra/reference/arith-generic.md)
for arithmetical operations. Use
[`roll`](https://rspatial.github.io/terra/reference/roll.md) for rolling
functions.

## Arguments

- x:

  SpatRaster

- base:

  a positive or complex number: the base with respect to which
  logarithms are computed

- digits:

  Number of digits for rounding

- fun:

  character. Math function name

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster or SpatExtent

## Examples

``` r
r1 <- rast(ncols=10, nrows=10)
v <- runif(ncell(r1))
v[10:20] <- NA
values(r1) <- v
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)
r <- c(r1, r2)

s <- sqrt(r)
# same as 
math(r, "sqrt")
#> class       : SpatRaster 
#> size        : 10, 10, 2  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       :     lyr.1, lyr.1 
#> min values  : 0.1272139,   0.1 
#> max values  : 0.9963790,   1.0 

round(s, 1)
#> class       : SpatRaster 
#> size        : 10, 10, 2  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.1 
#> min values  :   0.1,   0.1 
#> max values  :   1.0,   1.0 

cumsum(r)
#> class       : SpatRaster 
#> size        : 10, 10, 2  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       :      lyr.1,      lyr.1 
#> min values  : 0.01618339, 0.02618339 
#> max values  : 0.99277109, 1.84138276 
```
