# Disaggregate raster cells or vector geometries

`SpatRaster`: Create a SpatRaster with a higher resolution (smaller
cells). The values in the new SpatRaster are the same as in the larger
original cells.

`SpatVector`: Separate multi-objects (points, lines, polygons) into
single objects; or further into segments (for lines or polygons).

## Usage

``` r
# S4 method for class 'SpatRaster'
disagg(x, fact, method="near", filename="", ...)

# S4 method for class 'SpatVector'
disagg(x, segments=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- fact:

  positive integer. Aggregation factor expressed as number of cells in
  each direction (horizontally and vertically). Or two integers
  (horizontal and vertical aggregation factor) or three integers (when
  also aggregating over layers)

- method:

  character. Either "near" for nearest or "bilinear" for bilinear
  interpolation

- segments:

  logical. Should (poly-)lines or polygons be disaggregated into their
  line-segments?

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`aggregate`](https://rspatial.github.io/terra/reference/aggregate.md),
[`resample`](https://rspatial.github.io/terra/reference/resample.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(ncols=10, nrows=10)
rd <- disagg(r, fact=c(10, 2))
ncol(rd)
#> [1] 20
nrow(rd)
#> [1] 100
values(r) <- 1:ncell(r)
rd <- disagg(r, fact=c(4, 2))
```
