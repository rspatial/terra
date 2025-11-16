# Transpose

Transpose a SpatRaster or SpatVector

## Usage

``` r
# S4 method for class 'SpatRaster'
t(x)

# S4 method for class 'SpatVector'
t(x)

# S4 method for class 'SpatRaster'
trans(x, filename="", ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`flip`](https://rspatial.github.io/terra/reference/flip.md)`, `[`rotate`](https://rspatial.github.io/terra/reference/rotate.md)

## Examples

``` r
r <- rast(nrows=18, ncols=36)
values(r) <- 1:ncell(r)
tr1 <- t(r)
tr2 <- trans(r)
ttr <- trans(tr2)
```
