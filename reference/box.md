# draw a box

Similar to [`box`](https://rdrr.io/r/graphics/box.html) allowing adding
a box around a map. This function will place the box around the mapped
area.

## Usage

``` r
add_box(...)
```

## Arguments

- ...:

  arguments passed to
  [`lines`](https://rspatial.github.io/terra/reference/lines.md)

## See also

[`add_legend`](https://rspatial.github.io/terra/reference/legend.md),
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`add_mtext`](https://rspatial.github.io/terra/reference/add_mtext.md)

## Examples

``` r
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
add_box(col="red", lwd=3, xpd=TRUE)
```
