# rescale

Rescale a SpatVector or SpatRaster. This may be useful to make small
[`inset`](https://rspatial.github.io/terra/reference/inset.md) maps or
for georeferencing.

## Usage

``` r
# S4 method for class 'SpatRaster'
rescale(x, fx=0.5, fy=fx, x0, y0)

# S4 method for class 'SpatVector'
rescale(x, fx=0.5, fy=fx, x0, y0)
```

## Arguments

- x:

  SpatVector or SpatRaster

- fx:

  numeric \> 0. The horizontal scaling factor

- fy:

  numeric \> 0. The vertical scaling factor

- x0:

  numeric. x-coordinate of the center of rescaling. If missing, the
  center of the extent of `x` is used

- y0:

  numeric. y-coordinate of the center of rescaling. If missing, the
  center of the extent of `x` is used

## Value

Same as `x`

## See also

[`t`](https://rspatial.github.io/terra/reference/transpose.md),
[`shift`](https://rspatial.github.io/terra/reference/shift.md),
[`flip`](https://rspatial.github.io/terra/reference/flip.md),
[`rotate`](https://rspatial.github.io/terra/reference/rotate.md),
[`inset`](https://rspatial.github.io/terra/reference/inset.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- rescale(v, 0.2)
plot(v)
lines(w, col="red")
```
