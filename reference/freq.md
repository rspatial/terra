# Frequency table

Frequency table of the values of a SpatRaster. `NA`s are not counted
unless `value=NA`.

You can provide a SpatVector or additional SpatRaster to define zones
for which to do tabulations.

## Usage

``` r
# S4 method for class 'SpatRaster'
freq(x, digits=0, value=NULL, bylayer=TRUE, usenames=FALSE, 
    zones=NULL, wide=FALSE, touches=FALSE)
```

## Arguments

- x:

  SpatRaster

- digits:

  integer. Used for rounding the values before tabulation. Ignored if
  `NA`

- value:

  numeric. An optional single value to only count the number of cells
  with that value. This value can be `NA`

- bylayer:

  logical. If `TRUE` tabulation is done by layer

- usenames:

  logical. If `TRUE` layers are identified by their names instead of
  their numbers Only relevant if `bylayer` is `TRUE`

- zones:

  SpatRaster or SpatVector to define zones for which the tabulation
  should be done

- wide:

  logical. Should the results by "wide" instead of "long"?

- touches:

  logical. If `TRUE`, all cells touched by lines or polygons will be
  included, not just those on the line render path, or whose center
  point is within the polygon. Only relevant if `zones` is a SpatVector

## Value

A `data.frame` with 3 columns (layer, value, count) unless
`bylayer=FALSE` in which case a`data.frame` with two columns is returned
(value, count).

## Author

Andrew Gene Brown, Robert J. Hijmans

## Examples

``` r
r <- rast(nrows=10, ncols=10)
set.seed(2)
values(r) <- sample(5, ncell(r), replace=TRUE)

freq(r)
#>   layer value count
#> 1     1     1    27
#> 2     1     2    15
#> 3     1     3    17
#> 4     1     4    17
#> 5     1     5    24

x <- c(r, r/3)
freq(x, bylayer=FALSE)
#>   value count
#> 1     0    27
#> 2     1    76
#> 3     2    39
#> 4     3    17
#> 5     4    17
#> 6     5    24
freq(x)
#>   layer value count
#> 1     1     1    27
#> 2     1     2    15
#> 3     1     3    17
#> 4     1     4    17
#> 5     1     5    24
#> 6     2     0    27
#> 7     2     1    49
#> 8     2     2    24

freq(x, digits=1)
#>    layer value count
#> 1      1   1.0    27
#> 2      1   2.0    15
#> 3      1   3.0    17
#> 4      1   4.0    17
#> 5      1   5.0    24
#> 6      2   0.3    27
#> 7      2   0.7    15
#> 8      2   1.0    17
#> 9      2   1.3    17
#> 10     2   1.7    24
freq(x, digits=-1)
#>   layer value count
#> 1     1     0    76
#> 2     1    10    24
#> 3     2     0   100

freq(x, value=5)
#>   layer value count
#> 1     1     5    24
#> 2     2     5     0
```
