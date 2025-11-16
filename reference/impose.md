# Impose the geometry of a SpatRaster to those in a SpatRasterCollection.

Warp the members of a SpatRasterCollection to match the geometry of a
SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRasterCollection'
impose(x, y, filename="", ...)
```

## Arguments

- x:

  SpatRasterCollection

- y:

  SpatRaster

- filename:

  character. Output filename

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`resample`](https://rspatial.github.io/terra/reference/resample.md)
