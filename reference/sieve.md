# Sieve filter

Apply a sieve filter. That is, remove "noise", by changing small clumps
of cells with a value that is different from the surrounding cells, to
the value of the largest neighboring clump.

Note that the numerical input values are truncated to integers.

## Usage

``` r
# S4 method for class 'SpatRaster'
sieve(x, threshold, directions=8, filename="", ...)
```

## Arguments

- x:

  SpatRaster, single layer with integer or categorical values

- threshold:

  positive integer. Only clumps smaller than this threshold will be
  removed

- directions:

  numeric to indicate which cells are connected. Either `4` to only
  consider the horizontal and vertical neighbors ("rook"), or `8` to
  consider the vertical, horizontal and diagonal neighbors

- filename:

  character. Output filename

- ...:

  Options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`focal`](https://rspatial.github.io/terra/reference/focal.md)

## Examples

``` r
r <- rast(nrows=18, ncols=18, xmin=0, vals=0, crs="local")
r[2, 5] <- 1
r[5:8, 2:3] <- 2
r[7:12, 10:15] <- 3
r[15:16, 15:18] <- 4
freq(r, bylayer=FALSE)
#>   value count
#> 1     0   271
#> 2     1     1
#> 3     2     8
#> 4     3    36
#> 5     4     8

x <- sieve(r, 8)
y <- sieve(r, 9)
```
