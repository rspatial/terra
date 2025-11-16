# Geometry type of a SpatVector

Get the geometry type (points, lines, or polygons) of a SpatVector. See
[`datatype`](https://rspatial.github.io/terra/reference/datatype.md) for
the data types of the fields (attributes, variables) of a SpatVector.

## Usage

``` r
# S4 method for class 'SpatVector'
geomtype(x)

# S4 method for class 'SpatVector'
is.points(x)

# S4 method for class 'SpatVector'
is.lines(x)

# S4 method for class 'SpatVector'
is.polygons(x)
```

## Arguments

- x:

  SpatVector

## Value

character

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

geomtype(v)
#> [1] "polygons"
is.polygons(v)
#> [1] TRUE
is.lines(v)
#> [1] FALSE
is.points(v)
#> [1] FALSE

names(v)
#> [1] "ID_1"   "NAME_1" "ID_2"   "NAME_2" "AREA"   "POP"   
datatype(v)
#> [1] "double" "string" "double" "string" "double" "double"
```
