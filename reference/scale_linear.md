# Scale values linearly

Linear scaling of raster cell values between a specified minimum and
maximum value.

## Usage

``` r
# S4 method for class 'SpatRaster'
scale_linear(x, min=0, max=1, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- min:

  minimum value to scale to

- max:

  maximum value to scale to

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`scale`](https://rspatial.github.io/terra/reference/scale.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
s1 <- scale_linear(r)
s2 <- scale_linear(r, 1, 10)
```
