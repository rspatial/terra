# Which cells are TRUE?

This method returns a single layer SpatRaster with cell values
indicating the first layer in the input that is `TRUE`. All numbers that
are not zero (or `FALSE`), are considered to be `TRUE`.

## Usage

``` r
# S4 method for class 'SpatRaster'
which.lyr(x)
```

## Arguments

- x:

  SpatRaster

## Value

SpatRaster

## See also

[`isTRUE`](https://rdrr.io/r/base/Logic.html),
[`which`](https://rdrr.io/r/base/which.html), See
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md)
for `which.min` and `which.max`

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- which.lyr(s > 100)
```
