# clamp time series data

clamp time-series datat that are S shaped. The value in layers before
the minimum value in a cell can be set to that minimum value, and the
value in layers after the maximum value for a cell can be set to that
maximum value.

## Usage

``` r
# S4 method for class 'SpatRaster'
clamp_ts(x, min=FALSE, max=TRUE, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- min:

  logical. If `TRUE` the time-series is clamped to the minimum value

- max:

  logical. If `TRUE` the time-series is clamped to the maximum value

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`clamp`](https://rspatial.github.io/terra/reference/clamp.md),
[`cummin`](https://rdrr.io/r/base/cumsum.html),
[`cummax`](https://rdrr.io/r/base/cumsum.html)

## Examples

``` r
sigm <- function(x) { .8 / (1 + exp(-(x-10))) + runif(length(x))/4 }
r <- rast(ncols=10, nrows=10, nlyr=50)
s <- seq(5.2, 15,.2)
set.seed(1)
values(r) <- t(replicate(100, sigm(s)))

x <- clamp_ts(r, TRUE, TRUE) 

plot(unlist(r[4]))
lines(unlist(x[4]))

```
