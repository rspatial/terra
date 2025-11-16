# Spatial selection

Geometrically subset SpatRaster or SpatVector (to be done) by drawing on
a plot (map).

Note that for many installations this does to work well on the default
RStudio plotting device. To work around that, you can first run
`dev.new(noRStudioGD = TRUE)` which will create a separate window for
plotting, then use
[`plot()`](https://rspatial.github.io/terra/reference/plot.md) followed
by `sel()` and click on the map. It may also help to set your RStudio
"Tools/Global Options/Appearance/Zoom" to 100

## Usage

``` r
# S4 method for class 'SpatRaster'
sel(x, ...)

# S4 method for class 'SpatVector'
sel(x, use="rec", show=TRUE, col="cyan", draw=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- use:

  character indicating what to draw. One of "rec" (rectangle) or "pol"
  (polygon)

- show:

  logical. If `TRUE` the selected geometries are shown on the map

- col:

  color to be used for drawing if `draw=TRUE`

- draw:

  logical. If `TRUE` the area drawn to select geometries is shown on the
  map

- ...:

  additional graphics arguments for drawing the selected geometries

## See also

[`crop`](https://rspatial.github.io/terra/reference/crop.md) and
[`intersect`](https://rspatial.github.io/terra/reference/intersect.md)
to make an intersection and
[`click`](https://rspatial.github.io/terra/reference/click.md) and
[`text`](https://rspatial.github.io/terra/reference/text.md) to see cell
values or geometry attributes.

Use [`draw`](https://rspatial.github.io/terra/reference/draw.md) to draw
a SpatExtent of SpatVector that you want to keep.

## Value

SpatRaster or SpatVector

## Examples

``` r
if (FALSE) { # \dontrun{
# select a subset of a SpatRaster
r <- rast(nrows=10, ncols=10)
values(r) <- 1:ncell(r)
plot(r)
s <- sel(r) # now click on the map twice

# plot the selection on a new canvas:
x11()
plot(s)

# vector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v)
x <- sel(v) # now click on the map twice
x
} # }
```
