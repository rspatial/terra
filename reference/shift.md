# Shift

Shift a SpatRaster, SpatVector or SpatExtent to another location.

## Usage

``` r
# S4 method for class 'SpatRaster'
shift(x, dx=0, dy=0, filename="", ...)

# S4 method for class 'SpatVector'
shift(x, dx=0, dy=0)

# S4 method for class 'SpatExtent'
shift(x, dx=0, dy=0)
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- dx:

  numeric. The shift in horizontal direction

- dy:

  numeric. The shift in vertical direction

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

Same as `x`

## See also

[`flip`](https://rspatial.github.io/terra/reference/flip.md),
[`rotate`](https://rspatial.github.io/terra/reference/rotate.md)

## Examples

``` r
r <- rast(xmin=0, xmax=1, ymin=0, ymax=1)
r <- shift(r, dx=1, dy=-1)

e <- ext(r)
shift(e, 5, 5)
#> SpatExtent : 6, 7, 4, 5 (xmin, xmax, ymin, ymax)
```
