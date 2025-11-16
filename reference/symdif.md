# Symmetrical difference

Symmetrical difference of polygons

## Usage

``` r
# S4 method for class 'SpatVector,SpatVector'
symdif(x, y)
```

## Arguments

- x:

  SpatVector

- y:

  SpatVector

## Value

SpatVector

## See also

[`erase`](https://rspatial.github.io/terra/reference/erase.md)

## Examples

``` r
p <- vect(system.file("ex/lux.shp", package="terra"))
b <- as.polygons(ext(6, 6.4, 49.75, 50))
#sd <- symdif(p, b)
#plot(sd, col=rainbow(12))
```
