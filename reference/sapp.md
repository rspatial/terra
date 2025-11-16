# Apply a terra function that takes only a single layer and returns a SpatRaster to all layers of a SpatRaster

Apply to all layers of a SpatRaster a function that only takes a single
layer SpatRaster and returns a SpatRaster (these are rare). In most
cases you can also use `lapply` or `sapply` for this.

Or apply the same method to each sub-dataset (SpatRaster) in a
SpatRasterDataset

## Usage

``` r
# S4 method for class 'SpatRaster'
sapp(x, fun, ..., filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRasterDataset'
sapp(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster or SpatRasterDataset

- fun:

  if `x` is a `SpatRaster`: a function that takes a SpatRaster argument
  and can be applied to each layer of `x` (e.g.
  [`terrain`](https://rspatial.github.io/terra/reference/terrain.md). if
  `x` is a `SpatRasterDataset`: a function that is applied to all layers
  of the SpatRasters in `x` (e.g. `mean`

- ...:

  additional arguments to be passed to `fun`

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

` `[`lapp`](https://rspatial.github.io/terra/reference/lapp.md)`, `[`app`](https://rspatial.github.io/terra/reference/app.md)`, `[`tapp`](https://rspatial.github.io/terra/reference/tapp.md)`, `[`lapply`](https://rdrr.io/r/base/lapply.html)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra")) + 1  

#SpatRasterDataset
sd <- sds(s*2, s/2)
y <- sapp(sd, mean)
z <- sapp(sd, function(i) 2 * mean(i))
```
