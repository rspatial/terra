# add a custom legend

Wrapper around [`legend`](https://rdrr.io/r/graphics/legend.html) that
allows adding a custom legend to a map using a keyword such as "topleft"
or "bottomright". This function will place the legend in the locations
within the mapped area as delineated by the axes.

## Usage

``` r
add_legend(x, y, xpd=TRUE, ...)
```

## Arguments

- x:

  The keyword to be used to position the legend (or the x coordinate)

- y:

  The y coordinate to be used to position the legend (is x is also a
  coordinate)

- xpd:

  logical. If `TRUE`, the legend can be added outside the map area

- ...:

  arguments passed to [`legend`](https://rdrr.io/r/graphics/legend.html)

## See also

[`add_box`](https://rspatial.github.io/terra/reference/box.md),
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`add_mtext`](https://rspatial.github.io/terra/reference/add_mtext.md)

## Examples

``` r
v <- vect(system.file("ex/lux.shp", package="terra"))
plot(v)
points(centroids(v), col="red")
legend("topleft", legend = "centroids", pch = 20, xpd=NA, bg="white", col="red")
add_legend("topright", legend = "centroids", pch = 20, col="red")
```
