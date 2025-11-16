# replace cell values

Substitute(replace) cell values of a SpatRaster with a new value. See
[`classify`](https://rspatial.github.io/terra/reference/classify.md) for
more complex/flexible replacement.

## Usage

``` r
# S4 method for class 'SpatRaster'
subst(x, from, to, others=NULL, raw=FALSE, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- from:

  numeric value(s). Normally a vector of the same length as \`to\`. If
  `x` has multiple layers, it can also be a matrix of numeric value(s)
  where `nrow(x) == length(to)`. In that case the output has a single
  layer, with values based on the combination of the values of the input
  layers

- to:

  numeric value(s). Normally a vector of the same length as \`from\`. If
  `x` has a single layer, it can also be a matrix of numeric value(s)
  where `nrow(x) == length(from)`. In that case the output has multiple
  layers, one for each column in `to`

- others:

  numeric. If not `NULL` all values that are not matched are set to this
  value. Otherwise they retain their original value.

- raw:

  logical. If `TRUE`, the values in from and to are the raw cell values,
  not the categorical labels. Only relevant if `is.factor(x)`

- filename:

  character. Output filename

- ...:

  Additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`classify`](https://rspatial.github.io/terra/reference/classify.md),
[`clamp`](https://rspatial.github.io/terra/reference/clamp.md)

## Examples

``` r
r <- rast(ncols=5, nrows=5, xmin=0, xmax=1, ymin=0, ymax=1, crs="")
r <- init(r, 1:6)
x <- subst(r, 3, 7)
x <- subst(r, 2:3, NA)
x <- subst(x, NA, 10)

# multiple output layers
z <- subst(r, 2:3, cbind(20,30))

# multiple input layers
rr <- c(r, r+1, r+2)
m <- rbind(c(1:3), c(3:5))
zz <- subst(rr, m, c(100, 200))
```
