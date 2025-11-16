# Cell level regression

Run a regression model for each cell of a SpatRaster. The independent
variable can either be defined by a vector, or another SpatRaster to
make it spatially variable.

## Usage

``` r
# S4 method for class 'SpatRaster,numeric'
regress(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatRaster,SpatRaster'
regress(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)
```

## Arguments

- y:

  SpatRaster

- x:

  SpatRaster or numeric (of the same length as `nlyr(x)`

- formula:

  regression formula in the general form of `y ~ x`. You can add
  additional terms such as `I(x^2)`

- na.rm:

  logical. Remove NA values?

- cores:

  positive integer. If `cores > 1`, a 'parallel' package cluster with
  that many cores is created and used. You can also supply a cluster
  object.

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- regress(s, 1:nlyr(s))
```
