# Quantiles of spatial data

Compute quantiles for each cell across the layers of a SpatRaster.

You can use use
[`global`](https://rspatial.github.io/terra/reference/global.md)`(x, fun=quantile)`
to instead compute quantiles across cells for each layer.

You can also use this method to compute quantiles of the numeric
variables of a SpatVector.

## Usage

``` r
# S4 method for class 'SpatRaster'
quantile(x, probs=seq(0, 1, 0.25), na.rm=FALSE, filename="", ...) 

# S4 method for class 'SpatVector'
quantile(x, probs=seq(0, 1, 0.25), ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- probs:

  numeric vector of probabilities with values in \[0,1\]

- na.rm:

  logical. If `TRUE`, `NA`'s are removed from `x` before the quantiles
  are computed

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster with layers representing quantiles

## See also

[`app`](https://rspatial.github.io/terra/reference/app.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
rr <- c(r/2, r, r*2)
qr <- quantile(rr)
qr
#> class       : SpatRaster 
#> size        : 77, 101, 5  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> names       :    q0, q0.25, q0.5, q0.75,  q1 
#> min values  :   0.0,   0.0,    0,     0,   0 
#> max values  : 127.5, 127.5,  255,   510, 510 

if (FALSE) { # \dontrun{
# same but slower
qa <- app(rr, quantile)
} # }

#quantile by layer instead of by cell
qg <- global(r, quantile)
```
