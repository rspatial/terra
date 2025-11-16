# add a grid to a map made with terra

Adaptation of [`grid`](https://rdrr.io/r/graphics/grid.html) that allows
adding a grid to a map. This function will place the grid in the
locations within the mapped area as delineated by the axes.

If you set the tick marks yourself, you can use
[`add_abline`](https://rspatial.github.io/terra/reference/add_abline.md)
to create a grid:

Also see
[`graticule`](https://rspatial.github.io/terra/reference/graticule.md)

## Usage

``` r
add_grid(nx=NULL, ny=nx, col="lightgray", lty="dotted", lwd=1)
```

## Arguments

- nx, ny:

  number of cells of the grid in x and y direction. When NULL, as per
  default, the grid aligns with the tick marks on the corresponding
  default axis (i.e., tickmarks as computed by axTicks). When NA, no
  grid lines are drawn in the corresponding direction

- col:

  character or (integer) numeric; color of the grid lines

- lty:

  character or (integer) numeric; line type of the grid lines

- lwd:

  non-negative numeric giving line width of the grid lines

## See also

[`graticule`](https://rspatial.github.io/terra/reference/graticule.md),
[`add_abline`](https://rspatial.github.io/terra/reference/add_abline.md),
[`add_legend`](https://rspatial.github.io/terra/reference/legend.md),
[`add_box`](https://rspatial.github.io/terra/reference/box.md),
`add_grid`,
[`add_mtext`](https://rspatial.github.io/terra/reference/add_mtext.md)

## Examples

``` r
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
add_grid()
```
