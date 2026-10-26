# Add labels to a map

Plots labels, that is a textual (rather than color) representation of
values, on top an existing plot (map).

## Usage

``` r
# S4 method for class 'SpatRaster'
text(x, labels, digits=0, halo=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...)

# S4 method for class 'SpatVector'
text(x, labels, halo=FALSE, inside=FALSE, fit=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- labels:

  character. Optional. Vector of labels with `length(x)` or a variable
  name from `names(x)`

- digits:

  integer. How many digits should be used?

- halo:

  logical. If `TRUE` a "halo" is printed around the text

- hc:

  character. The halo color

- hw:

  numeric. The halo width

- inside:

  logical. if \\x) is a SpatVector of polygons, should the (x,y)
  coordinates of the label be inside the geometry

- fit:

  logical. if \\x) is a SpatVector of polygons, should an attempt be
  made assure that the labels are entirely inside their corresponding
  polygon

- jitter:

  numeric. The amount of random noise used to adjust label positions,
  possibly avoiding overlaps. See argument 'factor' in
  [`jitter`](https://rdrr.io/r/base/jitter.html)

- xpd:

  logical or `NA`. If `FALSE`, the text is clipped to the plot region,
  if `TRUE`, the text is clipped to the figure region, and if `NA`, the
  text is clipped to the device region

- ...:

  additional arguments to pass to graphics function
  [`text`](https://rdrr.io/r/graphics/text.html)

## See also

[`text`](https://rdrr.io/r/graphics/text.html)`, `[`plot`](https://rspatial.github.io/terra/reference/plot.md)`, `[`halo`](https://rspatial.github.io/terra/reference/halo.md)

## Value

matrix with the label locations (invisibly)

## Examples

``` r
r <- rast(nrows=4, ncols=4)
values(r) <- 1:ncell(r)

plot(r)
text(r)

set.seed(123)
text(r, jitter = 2, col = "red", halo = TRUE)


plot(r)
text(r, halo=TRUE, hc="blue", col="white", hw=0.2)


plot(r, col=rainbow(16))
xy <- text(r, col=c("black", "white"), vfont=c("sans serif", "bold"), cex=2)

head(xy)
#>         x    y
#> [1,] -135 67.5
#> [2,]  -45 67.5
#> [3,]   45 67.5
#> [4,]  135 67.5
#> [5,] -135 22.5
#> [6,]  -45 22.5
```
