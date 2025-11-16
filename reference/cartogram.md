# Cartogram

Make a cartogram, that is, a map where the area of polygons is made
proportional to another variable. This can be a good way to map raw
count data (e.g. votes).

## Usage

``` r
# S4 method for class 'SpatVector'
cartogram(x, var, type="nc")
```

## Arguments

- x:

  SpatVector

- var:

  character. A variable name in `x`

- type:

  character. Cartogram type, one of "nc" (non-contiguous) or "circles"
  (dorling)

## Value

SpatVector

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md),
[`rescale`](https://rspatial.github.io/terra/reference/rescale.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$value <- 1:12
pnc <- cartogram(v, "value", "nc")
pcirc <- cartogram(v, "value", "circles")
plot(v, col="light gray", border="gray")
lines(pnc, col="red", lwd=2)
lines(pcirc, col="blue", lwd=2)

```
