# Linear units of the coordinate reference system

Get the linear units of the coordinate reference system (crs) of a
SpatRaster or SpatVector expressed in m. The value returned is used
internally to transform area and perimeter measures to meters. The value
returned for longitude/latitude crs is zero.

## Usage

``` r
# S4 method for class 'SpatRaster'
linearUnits(x)

# S4 method for class 'SpatVector'
linearUnits(x)
```

## Arguments

- x:

  SpatRaster or SpatVector

## Value

numeric (meter)

## See also

[`crs`](https://rspatial.github.io/terra/reference/crs.md)

## Examples

``` r
x <- rast()
crs(x) <- ""
linearUnits(x)
#> [1] NaN

crs(x) <- "+proj=longlat +datum=WGS84"
linearUnits(x)
#> [1] 0

crs(x) <- "+proj=utm +zone=1 +units=cm"
linearUnits(x)
#> [1] 0.01

crs(x) <- "+proj=utm +zone=1 +units=km"
linearUnits(x)
#> [1] 1000

crs(x) <- "+proj=utm +zone=1 +units=us-ft"
linearUnits(x)
#> [1] 0.3048006
```
