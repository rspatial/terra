# head and tail of a SpatRaster or SpatVector

Show the head (first values) or tail (last values) of a SpatRaster or of
the attributes of a SpatVector.

## Usage

``` r
head(x, ...) 
tail(x, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- ...:

  additional arguments passed on to other methods

## Value

matrix (SpatRaster) or data.frame (SpatVector)

## See also

`show`, [`geom`](https://rspatial.github.io/terra/reference/geometry.md)

## Examples

``` r
r <- rast(nrows=25, ncols=25)
values(r) <- 1:ncell(r)
head(r)
#>   lyr.1
#> 1     1
#> 2     2
#> 3     3
#> 4     4
#> 5     5
#> 6     6
tail(r)
#>   lyr.1
#> 1   620
#> 2   621
#> 3   622
#> 4   623
#> 5   624
#> 6   625
```
