# mosaic SpatRasters

Combine adjacent and (partly) overlapping SpatRasters to form a single
new SpatRaster. Values in overlapping cells are averaged (by default) or
can be computed with another function.

With the "blend" function, smooth gradients in overlappoing zones are
created by weighting each raster's contribution by the distance from the
cell to the raster's nearest edge.

The SpatRasters must have the same origin and spatial resolution.

If rasters do not overlap, or if their values do not need to be
computed, you can use
[`merge`](https://rspatial.github.io/terra/reference/merge.md) instead.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
mosaic(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRasterCollection,missing'
mosaic(x, fun="mean", filename="", ...)
```

## Arguments

- x:

  SpatRaster or SpatRasterCollection

- y:

  SpatRaster if `x` is a SpatRaster, otherwise missing

- ...:

  additional SpatRasters

- fun:

  character. One of "blend", "mean", "median", "min", "max", "modal",
  "sum", "first", "last"

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`merge`](https://rspatial.github.io/terra/reference/merge.md)

## Examples

``` r
x <- rast(xmin=-110, xmax=-60, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-95, xmax=-45, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-80, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- mosaic(x, y, z)

m2 <- mosaic(z, y, x)

# with a SpatRasterCollection

spc <- sprc(list(x, y, z))

m <- mosaic(spc)

b <- mosaic(spc, fun="blend")
```
