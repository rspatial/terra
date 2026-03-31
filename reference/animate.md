# Animate a map

Animate (sequentially plot) the layers of a SpatRaster, or a
SpatVectorCollection, or the variables or geometries of a SpatVector, to
create a movie.

## Usage

``` r
# S4 method for class 'SpatRaster'
animate(x, pause=0.25, main, range=NULL, maxcell=50000, n=1, ...)

# S4 method for class 'SpatVector'
animate(x, pause=0.25, main="", n=1, vars=NULL, range=NULL, add=NULL, ...)

# S4 method for class 'SpatVectorCollection'
animate(x, pause=0.25, n=1, vars=NULL, range=NULL, ext=NULL, add=NULL, ...)
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
  range of all layers is used for rasters, or all variables for vectors
  if they are all numeric. If `NA` the range of each individual layer is
  used

- maxcell:

  positive integer. Maximum number of cells to use for the plot. If
  `maxcell < ncell(x)`, `spatSample(type="regular")` is used before
  plotting

- n:

  integer \> 0. Number of plotting loops

- vars:

  numeric or character to indicate the variables to animate. If this is
  NULL, the geometries are animated instead

- ext:

  `SpatExtent` object (or an object for which
  [`ext`](https://rspatial.github.io/terra/reference/ext.md) returns
  one) setting the spatial extent of the plot.

- add:

  logical. Add the geometries to the current plot? When looping over
  geometries: if `TRUE`, add all geometries to the current plot. If
  `NULL`, `add` is set to `FALSE` for the first geometry and `TRUE` for
  the remaining ones. This argument is ignored when `vars` is not `NULL`

- ...:

  additional arguments passed to
  [`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Value

None

## Author

A. Márcia Barbosa, Robert J. Hijmans

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
animate(r, n=1)




v <- vect(system.file("ex/lux.shp", package="terra"))
animate(v[1:3, ], n=1)

# animate(v, vars=names(v), pause=0.7)

s <- svc(as.lines(v[3:5,]), v, v[1:3,], as.points(v))
#animate(s, col="blue", alpha=0.3, pause=0.7)

# you can save an animation to file like this
# animation::saveGIF(terra::animate(v), "animation.gif")
```
