# Animate a map

Animate (sequentially plot) the layers of a SpatRaster, or the variables
or geometries of a SpatVector, to create a movie.

## Usage

``` r
# S4 method for class 'SpatRaster'
animate(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...)
# S4 method for class 'SpatVector'
animate(x, pause=0.25, main="", n=1, vars=NULL, range=NULL, add=NULL, ...)
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

- vars:

  numeric or character to indicate the variables to animate. If this is
  NULL, the geometries are animated instead

- range. Two numbers to fix the range of values in the legend between
  plots:

- add:

  logical. Add the geometries to the current map? When looping over
  variables: `NULL` is equivalent to `TRUE`. When looping over
  geometries: if `TRUE` , add all geometries to the current plot. If
  `NULL`, `add` is set to `FALSE` for the first geometry and `TRUE` for
  the remaining ones.

- ...:

  Additional arguments passed to
  [`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Value

None

## Author

Márcia Barbosa, Robert J. Hijmans

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
animate(s, n=1)




v <- vect(system.file("ex/lux.shp", package="terra"))
animate(v, n=2)

animate(v, n=1, vars=names(v))

```
