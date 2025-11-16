# Check for rotation

Check if a SpatRaster is "rotated" and needs to be rectified before it
can be used

See [`rectify`](https://rspatial.github.io/terra/reference/rectify.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
is.rotated(x)
```

## Arguments

- x:

  SpatRaster

## Value

logical. One value for each raster data \*source\*

## See also

[`rectify`](https://rspatial.github.io/terra/reference/rectify.md)`, `[`is.flipped`](https://rspatial.github.io/terra/reference/is.flipped.md)

## Examples

``` r
r <- rast(nrows=10, ncols=10, vals=1:100)
is.rotated(r)
#> [1] FALSE
```
