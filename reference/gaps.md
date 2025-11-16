# Find gaps between polygons

Get the gaps between polygons of a SpatVector

## Usage

``` r
# S4 method for class 'SpatVector'
gaps(x)
```

## Arguments

- x:

  SpatVector

## Value

SpatVector

## See also

[`sharedPaths`](https://rspatial.github.io/terra/reference/sharedPaths.md),
[`topology`](https://rspatial.github.io/terra/reference/topology.md),
and [`fillHoles`](https://rspatial.github.io/terra/reference/fill.md) to
get or remove polygon holes

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
h <- convHull(v[-12], "NAME_1")
g <- gaps(h)
```
