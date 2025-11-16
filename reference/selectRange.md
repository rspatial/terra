# Select the values of a range of layers, as specified by cell values in another SpatRaster

Use a single layer SpatRaster to select cell values from different
layers in a multi-layer SpatRaster. The values of the SpatRaster to
select layers (`y`) should be whole numbers between `1` and `nlyr(x)`
(values outside this range are ignored).

See [`rapp`](https://rspatial.github.io/terra/reference/rapp.md) for
applying a function to a range of variable size.

See [`extract`](https://rspatial.github.io/terra/reference/extract.md)
for extraction of values by cell, point, or otherwise.

## Usage

``` r
# S4 method for class 'SpatRaster'
selectRange(x, y, z=1, repint=0, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatRaster. Cell values must be positive integers. They indicate the
  first layer to select for each cell

- z:

  positive integer. The number of layers to select

- repint:

  integer \> 1 and \< nlyr(x) allowing for repeated selection at a fixed
  interval. For example, if `x` has 36 layers, and the value of a cell
  in `y`=2 and `repint` = 12, the values for layers 2, 14 and 26 are
  returned

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`rapp`](https://rspatial.github.io/terra/reference/rapp.md),
[`tapp`](https://rspatial.github.io/terra/reference/tapp.md),
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- 1
s <- c(r, r+2, r+5)
s <- c(s, s)
set.seed(1)
values(r) <- sample(3, ncell(r), replace=TRUE)
x <- selectRange(s, r)

x <- selectRange(s, r, 3)
```
