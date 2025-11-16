# Coerce to a "raster" object

Implementation of the generic
[`as.raster`](https://rdrr.io/r/grDevices/as.raster.html) function to
create a "raster" (small r) object. Such objects can be used for
plotting with the
[`rasterImage`](https://rdrr.io/r/graphics/rasterImage.html) function.
NOT TO BE CONFUSED with the Raster\* (big R) objects defined by the
'raster' package!

## Usage

``` r
# S4 method for class 'SpatRaster'
as.raster(x, maxcell=500000, col)
```

## Arguments

- x:

  SpatRaster

- maxcell:

  positive integer. Maximum number of cells to use for the plot

- col:

  vector of colors. The default is `map.pal("viridis", 100)`

## Value

'raster' object

## Examples

``` r
r <- rast(ncols=3, nrows=3)
values(r) <- 1:ncell(r)
as.raster(r)
#>      [,1]      [,2]      [,3]     
#> [1,] "#440154" "#472B7A" "#3A528B"
#> [2,] "#2C718E" "#218F8D" "#27AD80"
#> [3,] "#5BC863" "#ADDC2F" "#FDE725"
```
