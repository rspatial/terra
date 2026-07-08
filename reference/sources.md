# Data sources of a SpatRaster

Get the data sources of a SpatRaster or SpatVector or related object.
Sources are either files (or similar resources) or "", meaning that they
are in memory. You can use `hasValues` to check if in-memory layers
actually have cell values.

## Usage

``` r
# S4 method for class 'SpatRaster'
sources(x, nlyr=FALSE, bands=FALSE)

# S4 method for class 'SpatVector'
sources(x)

# S4 method for class 'SpatRaster'
hasValues(x)

# S4 method for class 'SpatRaster'
inMemory(x, bylayer=FALSE)
```

## Arguments

- x:

  SpatRaster, SpatRasterCollection, SpatVector or SpatVectorProxy

- nlyr:

  logical. If `TRUE` for each source, the number of layers is returned

- bands:

  logical. If `TRUE` for each source, the "bands" used, that is, the
  layer number in the source file, are returned

- bylayer:

  logical. If `TRUE` a value is returned for each layer instead of for
  each source

## Value

A vector of filenames, or `""` when there is no filename, if `nlyr` and
`bands` are both `FALSE`. Otherwise a `data.frame`

## See also

[`toMemory`](https://rspatial.github.io/terra/reference/toMemory.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
s <- rast(r)
values(s) <- 1:ncell(s)
rs <- c(r,r,s,r)
sources(rs)
#> [1] "/home/runner/work/_temp/Library/terra/ex/elev.tif"
#> [2] "/home/runner/work/_temp/Library/terra/ex/elev.tif"
#> [3] ""                                                 
#> [4] "/home/runner/work/_temp/Library/terra/ex/elev.tif"
hasValues(r)
#> [1] TRUE
x <- rast()
hasValues(x)
#> [1] FALSE
```
