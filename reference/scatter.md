# Scatterplot of two SpatRaster layers

Scatterplot of the values of two SpatRaster layers

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
plot(x, y, maxcell=100000, warn=TRUE, nc, nr, 
   maxnl=16, smooth=FALSE, gridded=FALSE, ncol=25, nrow=25, ...)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatRaster

- maxcell:

  positive integer. Maximum number of cells to use for the plot

- nc:

  positive integer. Optional. The number of columns to divide the
  plotting device in (when plotting multiple layers)

- nr:

  positive integer. Optional. The number of rows to divide the plotting
  device in (when plotting multiple layers)

- maxnl:

  positive integer. Maximum number of layers to plot (for multi-layer
  objects)

- smooth:

  logical. If `TRUE` show a smooth scatterplot (see
  [`smoothScatter`](https://rdrr.io/r/graphics/smoothScatter.html)

- gridded:

  logical. If `TRUE` the scatterplot is gridded (counts by cells)

- warn:

  boolean. Show a warning if a sample of the pixels is used (for
  scatterplot only)

- ncol:

  positive integer. Number of columns for gridding

- nrow:

  positive integer. Number of rows for gridding

- ...:

  additional graphical arguments

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
plot(s[[1]], s[[2]])

plot(s, sqrt(s[[3:1]]))
```
