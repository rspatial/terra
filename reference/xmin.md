# Get or set single values of an extent

Get or set single values of an extent. Values can be set for a
SpatExtent or SpatRaster, but not for a SpatVector)

## Usage

``` r
# S4 method for class 'SpatExtent'
xmin(x)

# S4 method for class 'SpatExtent'
xmax(x)

# S4 method for class 'SpatExtent'
ymin(x)

# S4 method for class 'SpatExtent'
ymax(x)

# S4 method for class 'SpatRaster'
xmin(x)

# S4 method for class 'SpatRaster'
xmax(x)

# S4 method for class 'SpatRaster'
ymin(x)

# S4 method for class 'SpatRaster'
ymax(x)

# S4 method for class 'SpatVector'
xmin(x)

# S4 method for class 'SpatVector'
xmax(x)

# S4 method for class 'SpatVector'
ymin(x)

# S4 method for class 'SpatVector'
ymax(x)

# S4 method for class 'SpatRaster,numeric'
xmin(x) <- value

# S4 method for class 'SpatRaster,numeric'
xmax(x) <- value

# S4 method for class 'SpatRaster,numeric'
ymin(x) <- value

# S4 method for class 'SpatRaster,numeric'
ymax(x) <- value
```

## Arguments

- x:

  SpatRaster, SpatExtent, or SpatVector

- value:

  numeric

## Value

SpatExtent or numeric coordinate

## See also

[`ext`](https://rspatial.github.io/terra/reference/ext.md)

## Examples

``` r
r <- rast()
ext(r)
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
ext(c(0, 20, 0, 20))
#> SpatExtent : 0, 20, 0, 20 (xmin, xmax, ymin, ymax)

xmin(r)
#> [1] -180
xmin(r) <- 0
xmin(r)
#> [1] 0
```
