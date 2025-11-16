# SpatRaster image method

Plot (make a map of) the values of a SpatRaster via
[`image`](https://rdrr.io/r/graphics/image.html). See
[`plot`](https://rspatial.github.io/terra/reference/plot.md) if you need
more fancy options such as a legend.

## Usage

``` r
# S4 method for class 'SpatRaster'
image(x, y=1, maxcell=500000, ...)
```

## Arguments

- x:

  SpatRaster

- y:

  positive integer indicating the layer to be plotted, or a character
  indicating the name of the layer

- maxcell:

  positive integer. Maximum number of cells to use for the plot

- ...:

  additional arguments as for
  `graphics::`[`image`](https://rdrr.io/r/graphics/image.html)

## See also

[`plot`](https://rspatial.github.io/terra/reference/plot.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
image(r)

image(r, col=rainbow(24))
```
