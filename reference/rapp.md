# Range-apply

Apply a function to a range of the layers of a SpatRaster that varies by
cell. The range is specified for each cell with one or two SpatRasters
(arguments `first` and `last`). For either `first` or `last` you can use
a single number instead.

You cannot use single numbers for both `first` and `last` because in
that case you could use
[`app`](https://rspatial.github.io/terra/reference/app.md) or
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md),
perhaps
[`subset`](https://rspatial.github.io/terra/reference/subset.md)ting the
layers of a SpatRaster.

See
[`selectRange`](https://rspatial.github.io/terra/reference/selectRange.md)
to create a new SpatRaster by extracting one or more values starting at
a cell-varying layer.

## Usage

``` r
# S4 method for class 'SpatRaster'
rapp(x, first, last, fun, ..., allyrs=FALSE, fill=NA, 
        clamp=FALSE, circular=FALSE, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- first:

  SpatRaster or positive integer between 1 and nlyr(x), indicating the
  first layer in the range of layers to be considered

- last:

  SpatRaster or positive integer between 1 and nlyr(x), indicating the
  last layer in the range to be considered

- fun:

  function to be applied

- ...:

  additional arguments passed to `fun`

- allyrs:

  logical. If `TRUE`, values for all layers are passed to `fun` but the
  values outside of the range are set to `fill`

- fill:

  numeric. The fill value for the values outside of the range, for when
  `allyrs=TRUE`

- clamp:

  logical. If `FALSE` and the specified range is outside `1:nlyr(x)` all
  cells are considered `NA`. Otherwise, the invalid part of the range is
  ignored

- circular:

  logical. If `TRUE` the values are considered circular, such as the
  days of the year. In that case, if first \> last the layers used are
  c(first:nlyr(x), 1:last). Otherwise, the range would be considered
  invalid and `NA` would be returned

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`selectRange`](https://rspatial.github.io/terra/reference/selectRange.md),
[`app`](https://rspatial.github.io/terra/reference/app.md),
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md),
[`lapp`](https://rspatial.github.io/terra/reference/lapp.md),
[`tapp`](https://rspatial.github.io/terra/reference/tapp.md)

## Examples

``` r
r <- rast(ncols=9, nrows=9)
values(r) <- 1:ncell(r)
s <- c(r, r, r, r, r, r)
s <- s * 1:6
s[1:2] <- NA
start <- end <- rast(r)
start[] <- 1:3
end[]   <- 4:6
a <- rapp(s, start, end, fun="mean")
b <- rapp(s, start, 2, fun="mean")

# cumsum from start to nlyr(x). return all layers
r <- rapp(s, start, nlyr(s), cumsum, allyrs=TRUE, fill=0)
# return only the final value
rr <- rapp(s, start, nlyr(s), function(i) max(cumsum(i)))
```
