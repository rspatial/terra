# Apply a function to subsets of layers of a SpatRaster

Apply a function to subsets of layers of a SpatRaster (similar to
[`tapply`](https://rdrr.io/r/base/tapply.html) and
[`aggregate`](https://rdrr.io/r/stats/aggregate.html)). The layers are
combined based on the `index`.

The number of layers in the output SpatRaster equals the number of
unique values in `index` times the number of values that the supplied
function returns for a single vector of numbers.

For example, if you have a SpatRaster with 6 layers, you can use
`index=c(1,1,1,2,2,2)` and `fun=sum`. This will return a SpatRaster with
two layers. The first layer is the sum of the first three layers in the
input SpatRaster, and the second layer is the sum of the last three
layers in the input SpatRaster. Indices are recycled such that
`index=c(1,2)` would also return a SpatRaster with two layers (one based
on the odd layers (1,3,5), the other based on the even layers (2,4,6)).

The index can also be one of the following values to group by time
period (if `x` has the appropriate
[`time`](https://rspatial.github.io/terra/reference/time.md) values):
"years", "months", "yearmonths", "dekads", "yeardekads", "weeks" (the
ISO 8601 week number, see Details), "yearweeks", "days", "doy" (day of
the year), "7days" (seven-day periods starting at Jan 1 of each year),
"10days", or "15days". It can also be a function that makes groups from
time values.

See [`app`](https://rspatial.github.io/terra/reference/app.md) or
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md)
if you want to use a more efficient function that returns multiple
layers based on **all** layers in the SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
tapp(x, index, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- index:

  factor or numeric (integer). Vector of length `nlyr(x)` (shorter
  vectors are recycled) grouping the input layers. It can also be one of
  the following values: "years", "months", "yearmonths", "days", "week"
  (ISO 8601 week number), or "doy" (day of the year)

- fun:

  function to be applied. The following functions have been
  re-implemented in C++ for speed: "sum", "mean", "median", "modal",
  "which", "which.min", "which.max", "min", "max", "prod", "any", "all",
  "sd", "std", "first". To use the base-R function for say, "min", you
  could use something like `fun = \(i) min(i)`

- ...:

  additional arguments passed to `fun`

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

"week" follows the ISO 8601 definition. Weeks start on Monday. If the
week containing 1 January has four or more days in the new year, then it
is considered week "01". Otherwise, it is the last week of the previous
year (week "52" or "53", and the next week is week 1.

## See also

[`app`](https://rspatial.github.io/terra/reference/app.md),
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
s <- c(r, r, r, r, r, r)
s <- s * 1:6
b1 <- tapp(s, index=c(1,1,1,2,2,2), fun=sum)
b1
#> class       : SpatRaster 
#> size        : 10, 10, 2  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       :  X1,   X2 
#> min values  :   6,   15 
#> max values  : 600, 1500 
b2 <- tapp(s, c(1,2,3,1,2,3), fun=sum)
b2
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       :  X1,  X2,  X3 
#> min values  :   5,   7,   9 
#> max values  : 500, 700, 900 
```
