# Lagged differences

Compute the difference between consecutive layers in a SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
diff(x, lag=1, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- lag:

  positive integer indicating which lag to use

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
d <- diff(s)
```
