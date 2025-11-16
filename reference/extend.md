# Extend

Enlarge the spatial extent of a SpatRaster. See
[`crop`](https://rspatial.github.io/terra/reference/crop.md) if you
(also) want to remove rows or columns.

Note that you can only enlarge SpatRasters with entire rows and columns.
Therefore, the extent of the output SpatRaster may not be exactly the
same as the requested. Depending on argument `snap` it may be a bit
smaller or larger.

You can also enlarge a SpatExtent with this method, or with an algebraic
notation (see examples)

## Usage

``` r
# S4 method for class 'SpatRaster'
extend(x, y, snap="near", fill=NA, filename="", overwrite=FALSE, ...) 

# S4 method for class 'SpatExtent'
extend(x, y)
```

## Arguments

- x:

  SpatRaster or SpatExtent

- y:

  If `x` is a SpatRaster, `y` should be a SpatExtent, or an object from
  which it can be extracted (such as SpatRaster and SpatVector objects).
  Alternatively, you can provide one, two or four non-negative integers
  indicating the number of rows and columns that need to be added at
  each side (a single positive integer when the number of rows and
  columns to be added is equal; or 2 number (columns, rows), or four
  (left column, right column, bottom row, top row). If `x` is a
  SpatExtent, `y` should likewise be a numeric vector of 1, 2, or 4
  elements

- snap:

  character. One of "near", "in", or "out". Used to align `y` to the
  geometry of `x`

- fill:

  numeric. The value used to for the new raster cells

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster or SpatExtent

## See also

[`trim`](https://rspatial.github.io/terra/reference/trim.md),
[`crop`](https://rspatial.github.io/terra/reference/crop.md),
[`merge`](https://rspatial.github.io/terra/reference/merge.md),
[`ext`](https://rspatial.github.io/terra/reference/ext.md),
[`resample`](https://rspatial.github.io/terra/reference/resample.md),
[`elongate`](https://rspatial.github.io/terra/reference/elongate.md)

## Examples

``` r
r <- rast(xmin=-150, xmax=-120, ymin=30, ymax=60, ncols=36, nrows=18)
values(r) <- 1:ncell(r)
e <- ext(-180, -100, 40, 70)
re <- extend(r, e)

# extend with a number of rows and columns (at each side)
re2 <- extend(r, c(2,10))

# SpatExtent
e <- ext(r)
e
#> SpatExtent : -150, -120, 30, 60 (xmin, xmax, ymin, ymax)
extend(e, 10)
#> SpatExtent : -160, -110, 20, 70 (xmin, xmax, ymin, ymax)
extend(e, c(10, -10, 0, 20))
#> SpatExtent : -160, -110, 30, 80 (xmin, xmax, ymin, ymax)


# add 10 columns / rows on all sides
e + 10
#> SpatExtent : -160, -110, 20, 70 (xmin, xmax, ymin, ymax)
# double extent
e * 2
#> SpatExtent : -165, -105, 15, 75 (xmin, xmax, ymin, ymax)
# increase extent by 25%
e * 1.25
#> SpatExtent : -153.75, -116.25, 26.25, 63.75 (xmin, xmax, ymin, ymax)
```
