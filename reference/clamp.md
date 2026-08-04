# Clamp values

Clamp values to a minimum and maximum value. That is, all values below a
lower threshold value and above the upper threshold value become either
`NA`, or, if `values=TRUE`, become the threshold value

## Usage

``` r
# S4 method for class 'SpatRaster'
clamp(x, lower=-Inf, upper=Inf, values=TRUE, filename="", ...)

# S4 method for class 'numeric'
clamp(x, lower=-Inf, upper=Inf, values=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or numeric vector

- lower:

  numeric with the lowest acceptable value (you can specify a different
  value for each layer). Or a SpatRaster that has a single layer or the
  same number of layers as `x`

- upper:

  numeric with the highest acceptable value (you can specify a different
  value for each layer). Or a SpatRaster that has a single layer or the
  same number of layers as `x`

- values:

  logical. If `FALSE` values outside the clamping range become `NA`, if
  `TRUE`, they get the extreme values

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster A SpatRaster when `x` is a SpatRaster. A numeric vector when
`x` is numeric.

[`classify`](https://rspatial.github.io/terra/reference/classify.md)`, `[`subst`](https://rspatial.github.io/terra/reference/subst.md)

r \<- rast(ncols=10, nrows=10) values(r) \<- 1:ncell(r) rc \<- clamp(r,
25, 75) rc

spatial
