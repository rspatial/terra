# Perimeter or length

This method returns the length of lines or the perimeter of polygons.

When the coordinate reference system is not longitude/latitude, you may
get more accurate results by first transforming the data to
longitude/latitude with
[`project`](https://rspatial.github.io/terra/reference/project.md)

## Usage

``` r
# S4 method for class 'SpatVector'
perim(x)
```

## Arguments

- x:

  SpatVector

## Value

numeric (m)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
perim(v)
#>  [1] 117100.12  93477.28  84502.45  44919.14  85032.61  74708.05  57991.42
#>  [8]  81203.93  74443.82  95564.74  80618.76  70810.67
```
