# Is a SpatRaster is flipped

Check if a SpatRaster is "flipped" vertically, and may need to be
adjusted with
[`flip`](https://rspatial.github.io/terra/reference/flip.md) before it
can be used.

## Usage

``` r
# S4 method for class 'SpatRaster'
is.flipped(x)
```

## Arguments

- x:

  SpatRaster

## Value

logical. One value for each raster data \*source\*

## See also

[`flip`](https://rspatial.github.io/terra/reference/flip.md)`, `[`is.rotated`](https://rspatial.github.io/terra/reference/is.rotated.md)

## Examples

``` r
r <- rast(nrows=10, ncols=10)
is.flipped(r)
#> [1] FALSE
```
