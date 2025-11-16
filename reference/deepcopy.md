# Deep copy

Make a deep copy of a SpatRaster or SpatVector. This is occasionally
useful when using an in-place replacement function that does not make
copy, such as
[`set.ext`](https://rspatial.github.io/terra/reference/inplace.md).

## Usage

``` r
# S4 method for class 'SpatRaster'
deepcopy(x)

# S4 method for class 'SpatVector'
deepcopy(x)
```

## Arguments

- x:

  SpatRaster or SpatVector

## Value

Same as `x`

## Examples

``` r
r <- rast(ncols=10, nrows=10, nl=3)
x <- r
y <- deepcopy(r)
ext(r)
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
set.ext(x, c(0,10,0,10))
ext(x)
#> SpatExtent : 0, 10, 0, 10 (xmin, xmax, ymin, ymax)
ext(r)
#> SpatExtent : 0, 10, 0, 10 (xmin, xmax, ymin, ymax)
ext(y)
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
```
