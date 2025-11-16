# Fill time gaps in a SpatRaster

Add empty layers in between existing layers such that the time step
between each layer is the same. See
[`approximate`](https://rspatial.github.io/terra/reference/approximate.md)
to estimate values for these layer (and other missing values)

## Usage

``` r
# S4 method for class 'SpatRaster'
fillTime(x, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- filename:

  character. Output filename

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`approximate`](https://rspatial.github.io/terra/reference/approximate.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
s <- c(r, r)
time(s) <- as.Date("2001-01-01") + c(0:2, 5:7)
time(s)
#> [1] "2001-01-01" "2001-01-02" "2001-01-03" "2001-01-06" "2001-01-07"
#> [6] "2001-01-08"
ss <- fillTime(s)
time(ss)
#> [1] "2001-01-01" "2001-01-02" "2001-01-03" "2001-01-04" "2001-01-05"
#> [6] "2001-01-06" "2001-01-07" "2001-01-08"

a <- approximate(ss)
```
