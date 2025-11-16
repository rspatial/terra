# add vertical and/or horizontal lines to a map made with terra

Adaptation of [`abline`](https://rdrr.io/r/graphics/abline.html) that
allows adding a horizonal or vertical lines to a map. This function will
place the lines in the locations within the mapped area as delineated by
the axes. It is meant to be used when you specifiy your own tick marks,
such that
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md) does
not work.

Also see
[`graticule`](https://rspatial.github.io/terra/reference/graticule.md)

## Usage

``` r
add_abline(h=NULL, v=NULL, ...)
```

## Arguments

- h:

  the y-value(s) for horizontal line(s)

- v:

  the x-value(s) for vertical line(s)

- ...:

  additional graphical parameters for drawing lines

## See also

[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`graticule`](https://rspatial.github.io/terra/reference/graticule.md),
[`add_legend`](https://rspatial.github.io/terra/reference/legend.md),
[`add_box`](https://rspatial.github.io/terra/reference/box.md),
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`add_mtext`](https://rspatial.github.io/terra/reference/add_mtext.md)

## Examples

``` r
v <- vect(system.file("ex/lux.shp", package="terra"))
atx <- seq(xmin(v), xmax(v), .1)
aty <- seq(ymin(v), ymax(v), .1)
plot(v, pax=list(xat=atx, yat=aty), ext=ext(v)+.2)
add_abline(h=aty, v=atx, lty=2, col="gray")
```
