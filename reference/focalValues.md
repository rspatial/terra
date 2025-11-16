# Get focal values

Get a matrix in which each row had the focal values of a cell. These are
the values of a cell and a rectangular window around it.

## Usage

``` r
# S4 method for class 'SpatRaster'
focalValues(x, w=3, row=1, nrows=nrow(x), fill=NA)
```

## Arguments

- x:

  SpatRaster or SpatVector

- w:

  window. The window can be defined as one (for a square) or two odd
  numbers (row, col); or with an odd sized matrix

- row:

  positive integer. Row number to start from, should be between 1 and
  nrow(x)

- nrows:

  positive integer. How many rows?

- fill:

  numeric used as values for imaginary cells outside the raster

## Value

matrix

## Examples

``` r
r <- rast(ncol=4, nrow=4, crs="+proj=utm +zone=1 +datum=WGS84")
values(r) <- 1:ncell(r)
focalValues(r)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]   NA   NA   NA   NA    1    2   NA    5    6
#>  [2,]   NA   NA   NA    1    2    3    5    6    7
#>  [3,]   NA   NA   NA    2    3    4    6    7    8
#>  [4,]   NA   NA   NA    3    4   NA    7    8   NA
#>  [5,]   NA    1    2   NA    5    6   NA    9   10
#>  [6,]    1    2    3    5    6    7    9   10   11
#>  [7,]    2    3    4    6    7    8   10   11   12
#>  [8,]    3    4   NA    7    8   NA   11   12   NA
#>  [9,]   NA    5    6   NA    9   10   NA   13   14
#> [10,]    5    6    7    9   10   11   13   14   15
#> [11,]    6    7    8   10   11   12   14   15   16
#> [12,]    7    8   NA   11   12   NA   15   16   NA
#> [13,]   NA    9   10   NA   13   14   NA   NA   NA
#> [14,]    9   10   11   13   14   15   NA   NA   NA
#> [15,]   10   11   12   14   15   16   NA   NA   NA
#> [16,]   11   12   NA   15   16   NA   NA   NA   NA
```
