# Rectify a SpatRaster

Rectify a rotated SpatRaster into a non-rotated object

## Usage

``` r
# S4 method for class 'SpatRaster'
rectify(x, method="bilinear", aoi=NULL, snap=TRUE,
         filename="", ...)
```

## Arguments

- x:

  SpatRaster to be rectified

- method:

  character. Method used to for resampling. See
  [`resample`](https://rspatial.github.io/terra/reference/resample.md)

- aoi:

  SpatExtent or SpatRaster to crop `x` to a smaller area of interest;
  Using a SpatRaster allowing to set the exact output extent and output
  resolution

- snap:

  logical. If `TRUE`, the origin and resolution of the output are the
  same as would the case when `aoi = NULL`. Only relevant if `aoi` is a
  `SpatExtent`

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`is.rotated`](https://rspatial.github.io/terra/reference/is.rotated.md)

## Value

SpatRaster
