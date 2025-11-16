# Align a SpatExtent

Align an SpatExtent with a SpatRaster This can be useful to create a new
SpatRaster with the same origin and resolution as an existing
SpatRaster. Do not use this to force data to match that really does not
match (use e.g.
[`resample`](https://rspatial.github.io/terra/reference/resample.md) or
(dis)aggregate for this).

It is also possible to align a SpatExtent to a clean divisor.

## Usage

``` r
# S4 method for class 'SpatExtent,SpatRaster'
align(x, y, snap="near")

# S4 method for class 'SpatExtent,numeric'
align(x, y)
```

## Arguments

- x:

  SpatExtent

- y:

  SpatRaster or numeric

- snap:

  Character. One of "near", "in", or "out", to determine in which
  direction the extent should be aligned. To the nearest border, inwards
  or outwards

## Value

SpatExtent

## See also

[`ext`](https://rspatial.github.io/terra/reference/ext.md),
[`draw`](https://rspatial.github.io/terra/reference/draw.md)

## Examples

``` r
r <- rast()
e <- ext(-10.1, 9.9, -20.1, 19.9)
ea <- align(e, r)
e
#> SpatExtent : -10.1, 9.9, -20.1, 19.9 (xmin, xmax, ymin, ymax)
ext(r)
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
ea
#> SpatExtent : -10, 10, -20, 20 (xmin, xmax, ymin, ymax)

align(e, 0.5)
#> SpatExtent : -10.5, 10, -20.5, 20 (xmin, xmax, ymin, ymax)
```
