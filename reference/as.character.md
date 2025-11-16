# Create a text representation of (the skeleton of) an object

Create a text representation of (the skeleton of) an object

## Usage

``` r
# S4 method for class 'SpatExtent'
as.character(x)

# S4 method for class 'SpatRaster'
as.character(x)
```

## Arguments

- x:

  SpatRaster

## Value

character

## Examples

``` r
r <- rast()
ext(r)
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
ext(c(0, 20, 0, 20))
#> SpatExtent : 0, 20, 0, 20 (xmin, xmax, ymin, ymax)
```
