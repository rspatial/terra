# Add labels to a map

Plots labels, that is a textual (rather than color) representation of
values, on top an existing plot (map).

## Usage

``` r
# S4 method for class 'SpatRaster'
text(x, labels, digits=0, halo=FALSE, hc="white", hw=0.1, jitter=0, ...)

# S4 method for class 'SpatVector'
text(x, labels, halo=FALSE, inside=FALSE, hc="white", hw=0.1, jitter=0, ...)
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

  logical. Should the text always be placed inside one the
  sub-geometries?

- jitter:

  numeric. The amount of random noise used to adjust label positions,
  possibly avoiding overlaps. See argument 'factor' in
  [`jitter`](https://rdrr.io/r/base/jitter.html)

- ...:

  additional arguments to pass to graphics function
  [`text`](https://rdrr.io/r/graphics/text.html)

## See also

[`text`](https://rdrr.io/r/graphics/text.html)`, `[`plot`](https://rspatial.github.io/terra/reference/plot.md)`, `[`halo`](https://rspatial.github.io/terra/reference/halo.md)

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
text(r, col=c("black", "white"), vfont=c("sans serif", "bold"), cex=2)
```
