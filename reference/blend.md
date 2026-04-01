# Distance-weighted raster blending

Combine rasters with smooth gradients in overlapping zones by weighting
each raster's contribution by the distance from the cell to the raster's
nearest edge. The result is independent of the order of the rasters.

The SpatRasters must have the same origin and spatial resolution.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
blend(x, y, ..., filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRasterCollection,missing'
blend(x, filename="", ...)
```

## Arguments

- x:

  SpatRaster or SpatRasterCollection

- y:

  SpatRaster if `x` is a SpatRaster, otherwise missing

- ...:

  additional SpatRasters if `x` is a SpatRaster, otherwise options for
  writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

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

[`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md),
[`merge`](https://rspatial.github.io/terra/reference/merge.md)

## Examples

``` r
x <- rast(xmin=-110, xmax=-60, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-95, xmax=-45, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-80, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- blend(x, y, z)
#> Error in blend(x, y, z): could not find function "blend"

# with many SpatRasters, make a SpatRasterCollection from a list
rlist <- list(x, y, z)
rsrc <- sprc(rlist)

m2 <- blend(rsrc)
#> Error in blend(rsrc): could not find function "blend"
```
