# force counter-clockwise polygons

Assure that the nodes of outer rings of polygons are in
counter-clockwise order.

## Usage

``` r
# S4 method for class 'SpatVector'
forceCCW(x)
```

## Arguments

- x:

  SpatVector of polygons

## Value

SpatVector

## Examples

``` r
p <- vect("POLYGON ((2 45, 2 55, 18 55, 18 45, 2 45))")
pcc <- forceCCW(p)
geom(pcc, wkt=TRUE)
#> [1] "POLYGON ((2 45, 18 45, 18 55, 2 55, 2 45))"
```
