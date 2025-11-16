# Interpolate points using a moving window

Interpolate points within a moving window using inverse distance
weighting. The maximum number of points used can be restricted,
optionally by selecting the nearest points.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatVector'
interpIDW(x, y, field, radius, power=2, smooth=0,
       maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...)

# S4 method for class 'SpatRaster,matrix'
interpIDW(x, y, radius, power=2, smooth=0, 
       maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatVector or matrix with three columns (x,y,z)

- field:

  character. field name in SpatVector `y`

- radius:

  numeric. The radius of the circle (single number). If `near=FALSE`, it
  is also possible to use two or three numbers. Two numbers are
  interpreted as the radii of an ellipse (x and y-axis). A third number
  should indicated the desired, counter clockwise, rotation of the
  ellipse (in degrees)

- power:

  numeric. Weighting power

- smooth:

  numeric. Smoothing parameter

- minPoints:

  numeric. The minimum number of points to use. If fewer points are
  found in a search ellipse it is considered empty and the fill value is
  returned

- maxPoints:

  numeric. The maximum number of points to consider in a search area.
  Additional points are ignored. If fewer points are found, the fill
  value is returned

- near:

  logical. Should the nearest points within the neighborhood be used if
  `maxPoints` is reached?

- fill:

  numeric. value to use to fill empty cells

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`rasterizeWin`](https://rspatial.github.io/terra/reference/rasterizeWin.md)`, `[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)`, `[`interpNear`](https://rspatial.github.io/terra/reference/interpNear.md)`, `[`interpolate`](https://rspatial.github.io/terra/reference/interpolate.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(ncol=100, nrow=100, crs="local", xmin=0, xmax=50, ymin=0, ymax=50)
set.seed(100)
x <- runif(25, 5, 45)
y <- runif(25, 5, 45)
z <- sample(25)
xyz <- cbind(x,y,z)

x <- interpIDW(r, xyz, radius=5, power=1, smooth=1, maxPoints=5)
```
