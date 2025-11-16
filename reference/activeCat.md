# Active category

Get or set the active category of a multi-categorical SpatRaster layer

## Usage

``` r
# S4 method for class 'SpatRaster'
activeCat(x, layer=1)
# S4 method for class 'SpatRaster'
activeCat(x, layer = 1) <- value
```

## Arguments

- x:

  SpatRaster

- layer:

  positive integer, the layer number or name

- value:

  positive integer or character, indicating which column in the
  categories to use. Note that when a number is used this index is zero
  based, and "1" refers to the second column. This is because the first
  column of the categories has the cell values, not categorical labels

## Value

integer

## See also

[`levels`](https://rspatial.github.io/terra/reference/factors.md),
[`cats`](https://rspatial.github.io/terra/reference/factors.md)

## Examples

``` r
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) + 10
d <- data.frame(id=11:13, cover=c("forest", "water", "urban"), letters=letters[1:3], value=10:12)
levels(r) <- d

activeCat(r)
#> [1] 1
activeCat(r) <- 3
activeCat(r)
#> [1] 3
```
