# Scale (gain) and offset

These functions can be used to get or set the scale (gain) and offset
parameters used to transform values when reading raster data from a
file. The parameters are applied to the raw values using the formula
below:

`value <- value * scale + offset`

The default value for scale is 1 and for offset is 0. 'scale' is
sometimes referred to as 'gain'.

Note that setting the scale and/or offset are intended to be used with
values that are stored in a file. When values are memory, assigning
scale or offset values will lead to the immediate computation of new
values; in such cases it would be clearer to use
[`Arith-methods`](https://rspatial.github.io/terra/reference/arith-generic.md).

## Usage

``` r
# S4 method for class 'SpatRaster'
scoff(x)

# S4 method for class 'SpatRaster'
scoff(x) <- value
```

## Arguments

- x:

  SpatRaster

- value:

  two-column matrix with scale (first column) and offset (second column)
  for each layer. Or `NULL` to remove all scale and offset values

## Value

matrix or changed SpatRaster

## Examples

``` r
r <- rast(system.file("ex/elev.tif", package="terra"))
minmax(r)
#>     elevation
#> min       141
#> max       547
scoff(r)
#>      scale offset
#> [1,]     1      0
r[4603]
#>   elevation
#> 1       279

scoff(r) <- cbind(10, 5)
minmax(r)
#>     elevation
#> min      1415
#> max      5475
scoff(r)
#>      scale offset
#> [1,]    10      5
r[4603]
#>   elevation
#> 1      2795
```
