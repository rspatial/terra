# Flip or reverse a raster

Flip the values of a SpatRaster by inverting the order of the rows
(`vertical=TRUE`) or the columns (`vertical=FALSE`).

`rev` is the same as a horizontal \*and\* a vertical flip.

## Usage

``` r
# S4 method for class 'SpatRaster'
flip(x, direction="vertical", filename="", ...)

# S4 method for class 'SpatVector'
flip(x, direction="vertical")

# S4 method for class 'SpatRaster'
rev(x)
```

## Arguments

- x:

  SpatRaster or SpatVector

- direction:

  character. Should (partially) match "vertical" to flip by rows, or
  "horizontal" to flip by columns

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`trans`](https://rspatial.github.io/terra/reference/transpose.md),
[`rotate`](https://rspatial.github.io/terra/reference/rotate.md)

## Examples

``` r
r <- rast(nrow=18, ncol=36)
m <- matrix(1:ncell(r), nrow=18)
values(r) <- as.vector(t(m))
rx <- flip(r, direction="h")

values(r) <- as.vector(m)
ry <- flip(r, direction="v")

v <- rev(r)
```
