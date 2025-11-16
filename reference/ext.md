# Create, get or set a SpatExtent

Get a SpatExtent of a SpatRaster, SpatVector, or other spatial objects.
Or create a SpatExtent from four numbers (xmin, xmax, ymin, ymax).

You can set the extent of a SpatRaster, but you cannot set the extent of
a SpatVector (see
[`rescale`](https://rspatial.github.io/terra/reference/rescale.md) for
that). See
[`set.ext`](https://rspatial.github.io/terra/reference/inplace.md) to
set the extent in place.

## Usage

``` r
# S4 method for class 'SpatRaster'
ext(x, cells=NULL)

# S4 method for class 'SpatVector'
ext(x)

# S4 method for class 'numeric'
ext(x, ..., xy=FALSE)

# S4 method for class 'SpatRaster,SpatExtent'
ext(x) <- value

# S4 method for class 'SpatRaster,numeric'
ext(x) <- value
```

## Arguments

- x:

  SpatRaster, SpatVector, a numeric vector of length four (xmin, xmax,
  ymin, ymax), a single numeric (xmin; see additional arguments
  under`...`), or missing (in which case the output is the global extent
  in lon-lat coordinates)

- cells:

  positive integer (cell) numbers to subset the extent to area covered
  by these cells

- value:

  SpatExtent, or numeric vector of length four (xmin, xmax, ymin, ymax)

- ...:

  if `x` is a single numeric value, additional numeric values for xmax,
  ymin, and ymax

- xy:

  logical. Set this to `TRUE` to indicate that coordinates are in (xmin,
  ymin, xmax, ymax) order, instead of in the terra standard order of
  (xmin, xmax, ymin, ymax)

## Value

A
[`SpatExtent`](https://rspatial.github.io/terra/reference/SpatExtent-class.md)
object.

## See also

[`xmin`](https://rspatial.github.io/terra/reference/xmin.md),
[`xmax`](https://rspatial.github.io/terra/reference/xmin.md),
[`ymin`](https://rspatial.github.io/terra/reference/xmin.md),
[`ymax`](https://rspatial.github.io/terra/reference/xmin.md)

## Examples

``` r
ext()
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)

r <- rast()
e <- ext(r)
as.vector(e)
#> xmin xmax ymin ymax 
#> -180  180  -90   90 
as.character(e)
#> [1] "ext(-180, 180, -90, 90)"

ext(r) <- c(0, 2.5, 0, 1.5)
r
#> class       : SpatRaster 
#> size        : 180, 360, 1  (nrow, ncol, nlyr)
#> resolution  : 0.006944444, 0.008333333  (x, y)
#> extent      : 0, 2.5, 0, 1.5  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
er <- ext(r)

round(er)
#> SpatExtent : 0, 3, 0, 2 (xmin, xmax, ymin, ymax)
# go "in"
floor(er)
#> SpatExtent : 0, 3, 0, 2 (xmin, xmax, ymin, ymax)
# go "out"
ceiling(er)
#> SpatExtent : 0, 2, 0, 1 (xmin, xmax, ymin, ymax)

ext(r) <- e
```
