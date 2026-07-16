# Visualize Directional Arrows on a SpatRaster

This function overlays directional arrows on a plotted
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.md)
object, typically used to represent vector fields such as wind
directions, drainage flow directions, or other angular data.

## Usage

``` r
arrows_on_rast(x="flowdir",x_base=NULL,unit=c("deg","rad","flowdir"),
                arrow_length=NA,clockwise=FALSE,
                angle_from_h_axis=0,
                angle_from_x_axis=angle_from_h_axis,
                length=0,...)
```

## Arguments

- x:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.md)
  object representing angles or directions (in degrees or radians or
  flowdirs).

- x_base:

  Optional. A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.md)
  object to use as a basemap if no plot is already present.

- unit:

  Character. Unit of the angle values in `x`. Either `"degrees"` or
  `"radians"` or `"flowdir"`. See
  [`terrain`](https://rspatial.github.io/terra/reference/terrain.md) for
  `flowdir` option.

- arrow_length:

  Numeric. Length of the arrows to be plotted.

- clockwise:

  Logical. Default is `FALSE`. If `TRUE`, angles measured in radians or
  degrees are considered clockwise.

- angle_from_x_axis:

  Anti-clockwise angle from x (or horizontal) axis to the reference
  direction. Default is 0.

- angle_from_h_axis:

  Anti-clockwise angle from x (or horizontal) axis to the reference
  direction. Default is 0, e.g., in the case of northing wind direction
  it is often used `angle_x_axis=90, unit="deg", clockwise=TRUE`.

- length,...:

  Additional graphical parameters passed to
  [`arrows`](https://rdrr.io/r/graphics/arrows.html).

## Value

Invisibly returns the coordinates and vector components used for
plotting.

## Author

Emanuele Cordano

## See also

[`arrows`](https://rdrr.io/r/graphics/arrows.html),
[`terrain`](https://rspatial.github.io/terra/reference/terrain.md)

## Examples

``` r

f <- system.file("ex/elev_vinschgau.tif", package="terra")
r <- rast(f)  |> aggregate(fact=2,fun=min)
d <- terrain(r, "flowdir")

plot(r,col=terrain.colors(10))
arrows_on_rast(d, unit="flowdir",col="black")

#> NULL
```
