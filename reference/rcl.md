# Combine row, column, and layer numbers

Get a matrix with the combination of row, column, and layer numbers

## Usage

``` r
# S4 method for class 'SpatRaster'
rcl(x, row=NULL, col=NULL, lyr=NULL)
```

## Arguments

- x:

  SpatRaster

- row:

  positive integer that are row number(s), a list thereof, or NULL for
  all rows

- col:

  as above for columns

- lyr:

  as above for layers

## Details

If a list is used for at least one of `row`, `col` or `lyr`, these are
evaluated in parallel. That is combinations are made for each list
element, not across list elements. If, in this case another argument is
not a list it has to have either length 1 (used for all cases) or have
the same length as the (longest) list, in which case the value is
coerced into a list with `as.list`

If multiple arguments are a list but they have different lengths,
theyare recycled to the longest list.

## Value

matrix

## See also

[`rowColCombine`](https://rspatial.github.io/terra/reference/xyCellFrom.md),
[`cellFromRowCol`](https://rspatial.github.io/terra/reference/xyCellFrom.md)

## Examples

``` r
x <- rast(ncol=5, nrow=5, nlyr=2)
values(x) <- 1:size(x)

rcl(x, 1, 2:3, 1:2)
#>      row col lyr
#> [1,]   1   2   1
#> [2,]   1   3   1
#> [3,]   1   2   2
#> [4,]   1   3   2

i <- rcl(x, 1, list(1:2, 3:4), 1:2)
i
#>      row col lyr
#> [1,]   1   1   1
#> [2,]   1   2   1
#> [3,]   1   3   2
#> [4,]   1   4   2

# get the values for these cells
x[i]
#> [1]  1  2 28 29
```
