# Animate a map

Animate (sequentially plot) the layers of a SpatRaster, or the
geometries of a SpatVector, to create a movie.

## Usage

``` r
# S4 method for class 'SpatRaster'
animate(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...)
# S4 method for class 'SpatVector'
animate(x, pause=0.25, main="", n=1, add=NULL, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- pause:

  numeric. How long should the pause be between layers?

- main:

  title for each layer. For SpatRaster, if not supplied, the z-value is
  used if available. Otherwise the names are used.

- range:

  numeric vector of length 2. Range of values to plot, If `NULL` the
  range of all layers is used. If `NA` the range of each individual
  layer is used

- maxcell:

  positive integer. Maximum number of cells to use for the plot. If
  `maxcell < ncell(x)`, `spatSample(type="regular")` is used before
  plotting

- n:

  integer \> 0. Number of plotting loops

- add:

  logical. If TRUE, add all geometries to the current plot. If left
  NULL, `add` is FALSE for the first geometry and TRUE for the remaining
  ones.

- ...:

  Additional arguments passed to
  [`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Value

None

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
animate(s, n=1)




v <- vect(system.file("ex/lux.shp", package="terra"))
animate(v, n=2)


animate(v[order(v$AREA), ])
```
