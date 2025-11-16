# Query by clicking on a map

Click on a map (plot) to get the coordinates or the values of a
SpatRaster or SpatVector at that location. For a SpatRaster you can also
get the coordinates and cell number of the location.

Note that for many installations this does to work well on the default
RStudio plotting device. To work around that, you can first run
`dev.new(noRStudioGD = TRUE)` which will create a separate window for
plotting, then use
[`plot()`](https://rspatial.github.io/terra/reference/plot.md) followed
by `click()` and click on the map. It may also help to set your RStudio
"Tools/Global Options/Appearance/Zoom" to 100

## Usage

``` r
# S4 method for class 'SpatRaster'
click(x, n=10, id=FALSE, xy=FALSE, cell=FALSE, type="p", show=TRUE, ...)

# S4 method for class 'SpatVector'
click(x,  n=10, id=FALSE, xy=FALSE, type="p", show=TRUE, ...)

# S4 method for class 'missing'
click(x, n=10, id=FALSE, type="p", show=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector, or missing

- n:

  number of clicks on the plot (map)

- id:

  logical. If `TRUE`, a numeric ID is shown on the map that corresponds
  to the row number of the output

- xy:

  logical. If `TRUE`, xy coordinates are included in the output

- cell:

  logical. If `TRUE`, cell numbers are included in the output

- type:

  one of "n", "p", "l" or "o". If "p" or "o" the points are plotted; if
  "l" or "o" they are joined by lines. See
  [`locator`](https://rdrr.io/r/graphics/locator.html)

- show:

  logical. Print the values after each click?

- ...:

  additional graphics parameters used if type != "n" for plotting the
  locations. See [`locator`](https://rdrr.io/r/graphics/locator.html)

## Value

The value(s) of `x` at the point(s) clicked on (or touched by the box
drawn). A `data.frame` with the value(s) of all layers of SpatRaster `x`
for the cell(s) clicked on; or with the attributes of the geometries of
SpatVector `x` that intersect with the box drawn).

## Note

The plot only provides the coordinates for a spatial query, the values
are read from the SpatRaster or SpatVector that is passed as an
argument. Thus, you can extract values from an object that has not been
plotted, as long as it spatially overlaps with the extent of the plot.

Unless the process is terminated prematurely values at most `n`
positions are determined. The identification process can be terminated,
depending on how you interact with R, by hitting Esc, or by clicking the
right mouse button and selecting "Stop" from the menu, or from the
"Stop" menu on the graphics window.

## See also

[draw](https://rspatial.github.io/terra/reference/draw.md)

## Examples

``` r
if (FALSE) { # \dontrun{
r <-rast(system.file("ex/elev.tif", package="terra"))
plot(r)
click(r, n=1)
## now click on the plot (map)
} # }
```
