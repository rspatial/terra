# Draw a polygon, line, extent, or points

Draw on a plot (map) to get a SpatVector or SpatExtent object for later
use. After calling the function, start clicking on the map. When you are
done, press `ESC`. You can also preset the maximum number of clicks.

Note that for many installations this does to work well on the default
RStudio plotting device. To work around that, you can first run
`dev.new(noRStudioGD = TRUE)` which will create a separate window for
plotting, then use
[`plot()`](https://rspatial.github.io/terra/reference/plot.md) followed
by `draw()` and clicking on the map. It may also help to set your
RStudio "Tools/Global Options/Appearance/Zoom" to 100

## Usage

``` r
# S4 method for class 'character'
draw(x="extent", col="red", lwd=2, id=FALSE, n=1000, xpd=TRUE, ...)
```

## Arguments

- x:

  character. The type of object to draw. One of "extent", "polygon",
  "line", or "points"

- col:

  the color to be used

- lwd:

  the width of the lines to be drawn

- id:

  logical. If `TRUE`, a numeric ID is shown on the map

- n:

  the maximum number of clicks (does not apply when `x=="extent"` in
  which case `n` is always 2)

- xpd:

  logical. If `TRUE`, you can draw outside the current plotting area

- ...:

  additional graphics arguments for drawing

## Value

SpatVector or SpatExtent

## See also

[`click`](https://rspatial.github.io/terra/reference/click.md)
