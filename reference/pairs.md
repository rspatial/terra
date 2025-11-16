# Pairs plot (matrix of scatterplots)

Pair plots of layers in a SpatRaster. This is a wrapper around graphics
function [`pairs`](https://rdrr.io/r/graphics/pairs.html).

## Usage

``` r
# S4 method for class 'SpatRaster'
pairs(x, hist=TRUE, cor=TRUE, use="pairwise.complete.obs", maxcells=100000, ...)
```

## Arguments

- x:

  SpatRaster

- hist:

  logical. If TRUE a histogram of the values is shown on the diagonal

- cor:

  logical. If TRUE the correlation coefficient is shown in the upper
  panels

- use:

  argument passed to the [`cor`](https://rdrr.io/r/stats/cor.html)
  function

- maxcells:

  integer. Number of pixels to sample from each layer of a large
  SpatRaster

- ...:

  additional arguments (graphical parameters)

## See also

[`boxplot`](https://rspatial.github.io/terra/reference/boxplot.md)`, `[`hist`](https://rspatial.github.io/terra/reference/hist.md)

## Examples

``` r
r <-rast(system.file("ex/elev.tif", package="terra"))
s <- c(r, 1/r, sqrt(r))
names(s) <- c("elevation", "inverse", "sqrt") 
pairs(s)


# to make indvidual histograms:
hist(r)

# or scatter plots:
plot(s[[1]], s[[2]])
```
