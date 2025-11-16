# Stretch

Linear or histogram equalization stretch of values in a SpatRaster.

For linear stretch, provide the desired output range (`minv` and `maxv`)
and the lower and upper bounds in the original data, either as quantiles
(`minq` and `maxq`, or as cell values (`smin` and `smax`). If `smin` and
`smax` are both not `NA`, `minq` and `maxq` are ignored.

For histogram equalization, these arguments are ignored, but you can
provide the desired scale of the output and the maximum number of cells
that is used to compute the histogram (empirical cumulative distribution
function).

## Usage

``` r
# S4 method for class 'SpatRaster'
stretch(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA,
    histeq=FALSE, scale=1, maxcell=500000, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- minv:

  numeric \>= 0 and smaller than maxv. lower bound of stretched value

- maxv:

  numeric \<= 255 and larger than maxv. upper bound of stretched value

- minq:

  numeric \>= 0 and smaller than maxq. lower quantile bound of original
  value. Ignored if smin is supplied

- maxq:

  numeric \<= 1 and larger than minq. upper quantile bound of original
  value. Ignored if smax is supplied

- smin:

  numeric \< smax. user supplied lower value for the layers, to be used
  instead of a quantile computed by the function itself

- smax:

  numeric \> smin. user supplied upper value for the layers, to be used
  instead of a quantile computed by the function itself

- histeq:

  logical. If `TRUE` histogram equalization is used instead of linear
  stretch

- scale:

  numeric. The scale (maximum value) of the output if `histeq=TRUE`

- maxcell:

  positive integer. The size of the regular sample used to compute the
  histogram

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(nc=10, nr=10)
values(r) <- rep(1:25, 4)
rs <- stretch(r)
s <- c(r, r*2)
sr <- stretch(s)
```
