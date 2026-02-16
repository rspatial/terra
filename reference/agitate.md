# Add noise to (jitter) a SpatVector of points

Add noise to (jitter) a SpatVector of points. Points are moved randomly
within a specified distance.

## Usage

``` r
# S4 method for class 'SpatVector'
agitate(x, maxdist, ...)
```

## Arguments

- x:

  SpatVector. All other objects are passed to
  [`jitter`](https://rdrr.io/r/base/jitter.html)

- maxdist:

  positive number. The maximum distance from the original location, in
  meter for lon/lat, otherwise in the units of the CRS

## Value

SpatVector

## Examples

``` r
## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- as.points(vect(f)[1])
p <- agitate(v, 2500)
```
