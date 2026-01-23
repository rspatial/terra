# Number of segments

Count the number of segments in a SpatVector of lines or polygons

## Usage

``` r
# S4 method for class 'SpatVector'
nseg(x)
```

## Arguments

- x:

  SpatVector

## Value

numeric

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
nseg(v)
#>  [1] 330 441 308 165 363 249 195 296 297 442 538 359
```
