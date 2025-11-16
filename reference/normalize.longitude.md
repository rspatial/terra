# normalize vector data that crosses the dateline

Normalize the longitude of geometries, move them if they are outside of
the -180 to 180 degrees range.

## Usage

``` r
# S4 method for class 'SpatVector'
normalize.longitude(x)
```

## Arguments

- x:

  SpatVector

## Value

SpatVector

## See also

[rotate](https://rspatial.github.io/terra/reference/rotate.md) for
SpatRaster

## Examples

``` r
p <- vect("POLYGON ((120 10, 230 75, 230 -75, 120 10))")
normalize.longitude(p)
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 1, 0  (geometries, attributes)
#>  extent      : -180, 180, -75, 75  (xmin, xmax, ymin, ymax)
#>  coord. ref. :  
```
