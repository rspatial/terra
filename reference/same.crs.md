# Compare coordinate reference systems

The function takes two coordinate reference system descriptions and
compares them for equality.

## Usage

``` r
same.crs(x, y)
```

## Arguments

- x:

  character, SpatRaster, SpatVector, CRS, or other object that returns
  something intelligible with`crs(x)`

- y:

  same types as for `x`

## Value

logical

## Examples

``` r
r <- rast()
same.crs(r, "+proj=longlat")
#> [1] TRUE

same.crs(r, "+proj=utm +zone=1")
#> [1] FALSE
```
