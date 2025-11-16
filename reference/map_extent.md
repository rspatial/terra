# Get the coordinates of the extent of a map

Helper function for creating custom map elements that are aligned with
the axes of a map (base plot created with a SpatRaster and/or
SpatVector). For example, you may need to know the coordinates for the
upper-left corner of a map to add some information there.

Unlike the standard base plot, terra keeps the axis aligned with the
data. For that reason you cannot use `par()$usr` to get these
coordinates.

The coordinates returned by this function are used in, for example,
[`add_legend`](https://rspatial.github.io/terra/reference/legend.md)
such that a legend can be automatically placed in the a particular
corner.

This function only returns meaningful results of the active plot
(canvas) was create with a call to `plot` with a SpatRaster or
SpatVector as first argument.

## Usage

``` r
map_extent()
```

## See also

[`add_legend`](https://rspatial.github.io/terra/reference/legend.md),
[`add_grid`](https://rspatial.github.io/terra/reference/grid.md),
[`add_box`](https://rspatial.github.io/terra/reference/box.md)

## Examples

``` r
r <- rast(xmin=0, xmax=10, ymin=0, ymax=10, res=1, vals=1:100)
plot(r)


map_extent()
#>  xmin xmax ymin ymax  geo
#>     0   10    0   10 TRUE
par()$usr
#> [1]  8.915711e-16  1.000000e+01 -5.484262e-01  1.054843e+01
```
