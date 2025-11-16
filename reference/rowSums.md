# row/col sums and means for SpatRaster

Sum or average values of SpatRaster layers by row or column.

## Usage

``` r
# S4 method for class 'SpatRaster'
rowSums(x, na.rm=FALSE, dims=1L, ...) 
# S4 method for class 'SpatRaster'
colSums(x, na.rm=FALSE, dims=1L, ...) 
# S4 method for class 'SpatRaster'
rowMeans(x, na.rm=FALSE, dims=1L, ...) 
# S4 method for class 'SpatRaster'
colMeans(x, na.rm=FALSE, dims=1L, ...)
```

## Arguments

- x:

  SpatRaster

- na.rm:

  logical. If `TRUE`, `NA` values are ignored

- dims:

  this argument is ignored

- ...:

  additional arguments (none implemented)

## Value

matrix

## See also

See [`global`](https://rspatial.github.io/terra/reference/global.md) for
summing all cells values

## Examples

``` r
r <- rast(ncols=2, nrows=5, nl=2, vals=1:20)
rowSums(r)
#>      lyr.1 lyr.2
#> [1,]     3    23
#> [2,]     7    27
#> [3,]    11    31
#> [4,]    15    35
#> [5,]    19    39
colSums(r)
#>      lyr.1 lyr.2
#> [1,]    25    75
#> [2,]    30    80
colMeans(r)
#>      lyr.1 lyr.2
#> [1,]     5    15
#> [2,]     6    16
```
