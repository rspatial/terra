# Apply a function to the cells of a SpatRaster

Apply a function to the values of each cell of a SpatRaster. Similar to
[`apply`](https://rdrr.io/r/base/apply.html) – think of each layer in a
SpatRaster as a column (or row) in a matrix.

This is generally used to summarize the values of multiple layers into
one layer; but this is not required.

`app` calls function `fun` with the raster data as first argument.
Depending on the function supplied, the raster data is represented as
either a matrix in which each layer is a column, or a vector
representing a cell. The function should return a vector or matrix that
is divisible by ncell(x). Thus, both "sum" and "rowSums" can be used,
but "colSums" cannot be used.

You can also apply a function `fun` across datasets by layer of a
`SpatRasterDataset`. In that case, summarization is by layer across
SpatRasters.

## Usage

``` r
# S4 method for class 'SpatRaster'
app(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRasterDataset'
app(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster or SpatRasterDataset

- fun:

  a function that operates on a vector or matrix. This can be a function
  that is defined in base-R or in a package, or a function you write
  yourself (see examples). Functions that return complex output (e.g. a
  list) may need to be wrapped in your own function to simplify the
  output to a vector or matrix. The following functions have been
  re-implemented in C++ for speed: "sum", "mean", "median", "modal",
  "which", "which.min", "which.max", "min", "max", "prod", "any", "all",
  "none", "sd", "std", "first". To use the base-R function for say,
  "min", you could use something like `fun=function(i) min(i)` or the
  equivalent `fun = \(i) min(i)`

- ...:

  additional arguments for `fun`. These are typically numerical
  constants. They should \*never\* be another SpatRaster

- cores:

  positive integer. If `cores > 1`, a 'parallel' package cluster with
  that many cores is created and used. You can also supply a cluster
  object. Ignored for functions that are implemented by terra in C++
  (see under fun)

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Details

To speed things up, parallelization is supported, but this is often not
helpful, and it may actually be slower. There is only a speed gain if
you have many cores (\> 8) and/or a very complex (slow) function `fun`.
If you write `fun` yourself, consider supplying a `cppFunction` made
with the Rcpp package instead (or go have a cup of tea while the
computer works for you).

## See also

[`lapp`](https://rspatial.github.io/terra/reference/lapp.md),
[`tapp`](https://rspatial.github.io/terra/reference/tapp.md),
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md),
[`roll`](https://rspatial.github.io/terra/reference/roll.md);
[`global`](https://rspatial.github.io/terra/reference/global.md) to
summarize the values of a single SpatRaster

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r, sqrt(r), r+50)
s <- app(x, fun=sum)
s
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : sum 
#> min value   :  53 
#> max value   : 260 
# for a few generic functions like 
# "sum", "mean", and "max" you can also do
sum(x)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : sum 
#> min value   :  53 
#> max value   : 260 

## SpatRasterDataset
sd <- sds(x, x*2, x/3)
a <- app(sd, max)
a
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.1, lyr.1 
#> min values  :     2,     2,   102 
#> max values  :   200,    20,   300 
# same as 
max(x, x*2, x/3)
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.1, lyr.1 
#> min values  :     2,     2,   102 
#> max values  :   200,    20,   300 
# and as (but slower)
b <- app(sd, function(i) max(i))


## also works for a single layer
f <- function(i) (i+1) * 2 * i + sqrt(i)
s <- app(r, f)
# same as above, but that is not memory-safe
# and has no filename argument 
s <- f(r)

if (FALSE) { # \dontrun{
#### multiple cores 
test0 <- app(x, sqrt) 
test1 <- app(x, sqrt, cores=2)

testfun <- function(i) { 2 * sqrt(i) }
test2 <- app(x, fun=testfun, cores =2)

## this fails because testfun is not exported to the nodes
# test3 <- app(x, fun=function(i) testfun(i), cores=2)
## to export it, add it as argument to fun
test3 <- app(x, fun=function(i, ff) ff(i), cores =3, ff=testfun)
} # }
```
