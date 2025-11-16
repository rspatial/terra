# modal value

Compute the mode for each cell across the layers of a SpatRaster. The
mode, or modal value, is the most frequent value in a set of values.

## Usage

``` r
# S4 method for class 'SpatRaster'
modal(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- ...:

  additional argument of the same type as `x` or numeric

- ties:

  character. Indicates how to treat ties. Either "random", "lowest",
  "highest", "first", or "NA"

- na.rm:

  logical. If `TRUE`, `NA` values are ignored. If `FALSE`, `NA` is
  returned if `x` has any `NA` values

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
r <- c(r/2, r, r*2)
m <- modal(r)
```
