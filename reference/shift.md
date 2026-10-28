# Shift

Shift a SpatRaster, SpatVector or SpatExtent to another location.

## Usage

``` r
# S4 method for class 'SpatRaster'
shift(x, dx=0, dy=0, filename="", to, ...)

# S4 method for class 'SpatVector'
shift(x, dx=0, dy=0, to)

# S4 method for class 'SpatExtent'
shift(x, dx=0, dy=0, to)
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- dx:

  numeric. The shift in horizontal direction

- dy:

  numeric. The shift in vertical direction

- to:

  SpatExtent, or an object from which an extent can be extracted (see
  [`ext`](https://rspatial.github.io/terra/reference/ext.md)); or two
  coordinates (x, y) for the new center of `x`. One of the two
  coordinates may be `NA` to shift only in the other direction. If `to`
  is used, `dx` and `dy` cannot be used

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

# center on another object or on coordinates
e2 <- ext(10, 12, 4, 8)
shift(e, to=e2)
#> SpatExtent : 10.5, 11.5, 5.5, 6.5 (xmin, xmax, ymin, ymax)
shift(e, to=c(11, 6))
#> SpatExtent : 10.5, 11.5, 5.5, 6.5 (xmin, xmax, ymin, ymax)
shift(e, to=c(11, NA))
#> SpatExtent : 10.5, 11.5, -1, 3.8163916471489756e-17 (xmin, xmax, ymin, ymax)
```
