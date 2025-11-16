# Shared paths

Get shared paths of line or polygon geometries. This can for geometries
in a single SpatVector, or between two SpatVectors

## Usage

``` r
# S4 method for class 'SpatVector'
sharedPaths(x, y=NULL)
```

## Arguments

- x:

  SpatVector of lines or polygons

- y:

  missing or SpatVector of lines or polygons

## Value

SpatVector

## See also

[`gaps`](https://rspatial.github.io/terra/reference/gaps.md),
[`topology`](https://rspatial.github.io/terra/reference/topology.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v, col="light gray")
text(v, halo=TRUE)

x <- sharedPaths(v)
lines(x, col="red", lwd=2)
text(x, col="blue", halo=TRUE, cex=0.8)

head(x)
#>   id1 id2
#> 1   1   2
#> 2   1   4
#> 3   1   5
#> 4   2   3
#> 5   2   4
#> 6   2   5

z <- sharedPaths(v[3,], v[12,])
```
