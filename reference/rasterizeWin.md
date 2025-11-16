# Rasterize points with a moving window

Rasterize points using a circle (or ellipse) as moving window. For each
raster cell, the points (`x, y`) that fall within the window centered on
that cell are considered. A function is used to compute a summary value
(e.g. "mean") for the values (`z`) associated with these points.

This can result in much smoother results compared to the standard
[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md)
method.

## Usage

``` r
# S4 method for class 'SpatVector,SpatRaster'
rasterizeWin(x, y, field, win="circle", pars, fun, ..., cvars=FALSE, 
        minPoints=1, fill=NA, filename="", wopt=list())

# S4 method for class 'data.frame,SpatRaster'
rasterizeWin(x, y, win="circle", pars, fun, ..., cvars=FALSE, 
          minPoints=1, fill=NA, filename="", wopt=list())
```

## Arguments

- x:

  SpatVector or matrix with at least three columns ((x, y) coordinates
  and a variable to be rasterized)

- y:

  SpatRaster

- field:

  character. field name in SpatVector `x` with the values to rasterize

- win:

  character to choose the window type. Can be "circle", "ellipse",
  "rectangle", or "buffer"

- pars:

  parameters to define the window. If `win="circle"` or `win="buffer"`,
  a single number to set the radius of the circle or the width of the
  buffer. If `win="ellipse"`, either two numbers (the x and y-axis) or
  three numbers the axes and a rotation (in degrees). If
  `win="rectangle"`, either two (width, height) or three (width, height)
  and the rotation in degrees. The unit of the radius/width/height/axis
  parameters is that of the coordinate reference system (it is not
  expressed as cells). That is, if you have a lon/lat crs, there is no
  conversion of degrees to meters or vice-versa.

- fun:

  function to summarize the values for each cell. If `cvars=FALSE`,
  functions must take a numeric vector and return (in all cases) one or
  more numbers. If `cvars=TRUE`, and multiple variables are used, the
  function must take a single argument (a data.frame with the names
  variables). For `win="circle"` and `win="ellipse"` there are two
  additional character values that can be used: `"distto"` (average
  distance to the points from the center of the cell) and
  `"distbetween"` (average distance between the points inside the
  window)

- ...:

  additional named arguments passed to `fun`

- minPoints:

  numeric. The minimum number of points to use. If fewer points are
  found in a search ellipse it is considered empty and the fill value is
  returned

- fill:

  numeric. value to use to fill cells with empty search areas

- cvars:

  logical. When using multiple fields, should `fun` operate on all of
  them at once? If not, `fun` is applied to each variable separately

- filename:

  character. Output filename

- wopt:

  list with additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`rasterize`](https://rspatial.github.io/terra/reference/rasterize.md),
[`rasterizeGeom`](https://rspatial.github.io/terra/reference/rasterizeGeom.md),
[`interpNear`](https://rspatial.github.io/terra/reference/interpNear.md),
[`interpIDW`](https://rspatial.github.io/terra/reference/interpIDW.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(ncol=100, nrow=100, crs="local", xmin=0, xmax=50, ymin=0, ymax=50)
set.seed(100)
x <- runif(50, 5, 45)
y <- runif(50, 5, 45)
z <- sample(50)
xyz <- data.frame(x,y,z)

r <- rasterizeWin(xyz, r, fun="count", pars=5)

rfuns <- c("count", "min", "max", "mean")
x <- lapply(rfuns, function(f) rasterizeWin(xyz, r, fun=f, pars=5))
names(x) <- rfuns 
x <- rast(x)
#plot(x)
```
