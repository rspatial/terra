# Weighted mean of layers

Compute the weighted mean for each cell of the layers of a SpatRaster.
The weights can be spatially variable or not.

## Usage

``` r
# S4 method for class 'SpatRaster,numeric'
weighted.mean(x, w, na.rm=FALSE, filename="", ...)

# S4 method for class 'SpatRaster,SpatRaster'
weighted.mean(x, w, na.rm=FALSE, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- w:

  A vector of weights (one number for each layer), or for spatially
  variable weights, a SpatRaster with weights (should have the same
  extent, resolution and number of layers as x)

- na.rm:

  Logical. Should missing values be removed?

- filename:

  character. Output filename

- ...:

  options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md),
[`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html)

## Examples

``` r
b <- rast(system.file("ex/logo.tif", package="terra"))   

# give least weight to first layer, most to last layer
wm1 <- weighted.mean(b, w=1:3)

# spatially varying weights
# weigh by column number
w1 <- init(b, "col")

# weigh by row number
w2 <- init(b, "row")
w <- c(w1, w2, w2)

wm2 <- weighted.mean(b, w=w)
```
