# spin a SpatVector

Spin (rotate) the geometry of a SpatVector.

## Usage

``` r
# S4 method for class 'SpatVector'
spin(x, angle, x0, y0)
```

## Arguments

- x:

  SpatVector

- angle:

  numeric. Angle of rotation in degrees

- x0:

  numeric. x-coordinate of the center of rotation. If missing, the
  center of the extent of `x` is used

- y0:

  numeric. y-coordinate of the center of rotation. If missing, the
  center of the extent of `x` is used

## Value

SpatVector

## See also

[`rescale`](https://rspatial.github.io/terra/reference/rescale.md),
[`t`](https://rspatial.github.io/terra/reference/transpose.md),
[`shift`](https://rspatial.github.io/terra/reference/shift.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- spin(v, 180)
plot(v)
lines(w, col="red")


# lower-right corner as center
e <- as.vector(ext(v))
x <- spin(v, 45, e[1], e[3])
```
