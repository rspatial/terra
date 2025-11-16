# units of SpatRaster or SpatRasterDataSet

Get or set the units of the layers of a SpatRaster or the datasets in a
SpatRasterDataSet.

## Usage

``` r
# S4 method for class 'SpatRaster'
units(x)

# S4 method for class 'SpatRaster'
units(x) <- value

# S4 method for class 'SpatRasterDataset'
units(x)

# S4 method for class 'SpatRasterDataset'
units(x) <- value
```

## Arguments

- x:

  SpatRaster

- value:

  character

## Value

character

## See also

[`time`](https://rspatial.github.io/terra/reference/time.md)`, `[`names`](https://rspatial.github.io/terra/reference/names.md)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   

units(s) <- c("m/s", "kg", "ha")
units(s)
#> [1] "m/s" "kg"  "ha" 
s
#> class       : SpatRaster 
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> colors RGB  : 1, 2, 3 
#> names       : red, green, blue 
#> min values  :   0,     0,    0 
#> max values  : 255,   255,  255 
#> unit        : m/s,    kg,   ha 

units(s) <- "kg"
units(s)
#> [1] "kg" "kg" "kg"
```
