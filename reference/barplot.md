# Bar plot of a SpatRaster

Create a barplot of the values of the first layer of a SpatRaster. For
large datasets a regular sample with a size of approximately `maxcells`
is used.

## Usage

``` r
# S4 method for class 'SpatRaster'
barplot(height, maxcell=1000000, digits=0, breaks=NULL, col, ...)
```

## Arguments

- height:

  SpatRaster

- maxcell:

  integer. To regularly subsample very large datasets

- digits:

  integer used to determine how to
  [`round`](https://rdrr.io/r/base/Round.html) the values before
  tabulating. Set to `NULL` or to a large number if you do not want any
  rounding

- breaks:

  breaks used to group the data as in
  [`cut`](https://rdrr.io/r/base/cut.html)

- col:

  a color generating function such as
  [`rainbow`](https://rdrr.io/r/grDevices/palettes.html) (the default),
  or a vector of colors

- ...:

  additional arguments for plotting as in
  [`barplot`](https://rdrr.io/r/graphics/barplot.html)

## See also

[`hist`](https://rspatial.github.io/terra/reference/hist.md)`, `[`boxplot`](https://rspatial.github.io/terra/reference/boxplot.md)

## Value

A numeric vector (or matrix, when `beside = TRUE`) of the coordinates of
the bar midpoints, useful for adding to the graph. See
[`barplot`](https://rdrr.io/r/graphics/barplot.html)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
barplot(r, digits=-1, las=2, ylab="Frequency")


op <- par(no.readonly = TRUE)
par(mai = c(1, 2, .5, .5))
barplot(r, breaks=10, col=c("red", "blue"), horiz=TRUE, digits=NULL, las=1)

par(op)
```
