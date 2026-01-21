# Cross-tabulate

Cross-tabulate the layers of a SpatRaster to create a contingency table.

## Usage

``` r
# S4 method for class 'SpatRaster,missing'
crosstab(x, digits=0, long=FALSE, useNA=FALSE)
```

## Arguments

- x:

  SpatRaster

- digits:

  integer. The number of digits for rounding the values before
  cross-tabulation

- long:

  logical. If `TRUE` the results are returned in 'long' format
  data.frame instead of a table

- useNA:

  logical, indicting if the table should includes counts of `NA` values

## Value

A table or data.frame

## See also

[`freq`](https://rspatial.github.io/terra/reference/freq.md),
[`zonal`](https://rspatial.github.io/terra/reference/zonal.md)

## Author

Andrew Gene Brown, Robert J. Hijmans

## Examples

``` r
r <- s <- rast(nc=5, nr=5)
set.seed(1)
values(r) <- runif(ncell(r)) * 2
values(s) <- runif(ncell(r)) * 3
x <- c(r, s)

crosstab(x)
#>      lyr.1.1
#> lyr.1 0 1 2 3
#>     0 1 1 4 0
#>     1 2 5 5 0
#>     2 0 2 4 1

rs <- r/s
r[1:5] <- NA
s[20:25] <- NA
x <- c(r, s, rs)
crosstab(x, useNA=TRUE, long=TRUE)
#>    lyr.1 lyr.1 lyr.1 n
#> 1      0     2     0 3
#> 2      0    NA     0 1
#> 3      0    NA     6 1
#> 4      1     0     4 1
#> 5      1     1     1 1
#> 6      1     1     2 1
#> 7      1     2     0 3
#> 8      1     2     1 1
#> 9      1    NA     0 1
#> 10     1    NA     1 1
#> 11     2     1     1 2
#> 12     2     2     1 2
#> 13     2    NA     1 2
#> 14    NA     0    19 1
#> 15    NA     1     0 2
#> 16    NA     1     1 1
#> 17    NA     3     1 1
```
