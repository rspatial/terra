# modal value

Compute the mode (most frequent value). For a SpatRaster this is done
for each cell across layers. For a vector (or other atomic object) it is
computed from the values of `x` and any additional arguments in `...`.

## Usage

``` r
# S4 method for class 'SpatRaster'
modal(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'ANY'
modal(x, ..., ties="random", na.rm=FALSE, freq=FALSE)
```

## Arguments

- x:

  SpatRaster or a vector (numeric, integer, logical, character, or
  factor)

- ...:

  additional argument of the same type as `x`, or numeric (SpatRaster
  method)

- ties:

  character. Indicates how to treat ties. Either "random", "lowest",
  "highest", "first", or "NA"

- na.rm:

  logical. If `TRUE`, `NA` values are ignored. If `FALSE`, `NA` is
  returned if `x` has any `NA` values

- freq:

  logical. If `TRUE`, the frequency of the modal value is returned
  instead of the value itself

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster, or a single value (or its frequency if `freq=TRUE`)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
r <- c(r/2, r, r*2)
m <- modal(r)

modal(c(1, 2, 2, 3, 1, 2), ties="lowest")
#> [1] 2
```
