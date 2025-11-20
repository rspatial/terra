# Transfer values of a SpatRaster to another one with a different geometry

resample transfers values between SpatRaster objects that do not align
(have a different origin and/or resolution). See
[`project`](https://rspatial.github.io/terra/reference/project.md) to
change the coordinate reference system (crs).

If the origin and extent of the input and output are the same, you
should consider using these other functions instead:
[`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md),
[`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md),
[`extend`](https://rspatial.github.io/terra/reference/extend.md) or
[`crop`](https://rspatial.github.io/terra/reference/crop.md).

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
resample(x, y, method, threads=FALSE, by_util=FALSE, filename="", ...)
```

## Arguments

- x:

  SpatRaster to be resampled

- y:

  SpatRaster with the geometry that `x` should be resampled to. You can
  also provide one or two positive numbers to set the resolution of the
  output raster relative to the input raster

- method:

  character. Method used for estimating the new cell values. One of:

  `bilinear`: bilinear interpolation (3x3 cell window). This is used by
  default if the first layer of `x` is not categorical

  `mean`: This can be a good choice with continuous variables if the
  output cells overlap with multiple input cells.

  `near`: nearest neighbor. This is used by default if the first layer
  of `x` is categorical. This method is not a good choice for continuous
  values.

  `modal`: The modal value. This can be a good choice for categorical
  rasters, if the output cells overlap with multiple input cells.

  `cubic`: cubic interpolation (5x5 cell window).

  `cubicspline`: cubic B-spline interpolation. (5x5 cell window).

  `lanczos`: Lanczos windowed sinc resampling. (7x7 cell window).

  `sum`: the weighted sum of all non-NA contributing grid cells.

  `min, q1, median, q3, max`: the minimum, first quartile, median, third
  quartile, or maximum value.

  `rms`: the root-mean-square value of all non-NA contributing grid
  cells.

- threads:

  logical. If `TRUE` multiple threads are used (faster for large files)

- by_util:

  logical. If `TRUE` the GDAL warp utility is used

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md),
[`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md),
[`crop`](https://rspatial.github.io/terra/reference/crop.md),
[`project`](https://rspatial.github.io/terra/reference/project.md)

## Examples

``` r
r <- rast(nrows=3, ncols=3, xmin=0, xmax=10, ymin=0, ymax=10)
values(r) <- 1:ncell(r)
s <- rast(nrows=25, ncols=30, xmin=1, xmax=11, ymin=-1, ymax=11)
x <- resample(r, s, method="bilinear")

opar <- par(no.readonly =TRUE)
par(mfrow=c(1,2))
plot(r)
plot(x)

par(opar)
```
