# Combine SpatRaster or SpatVector objects

With `c` you can:

– Combine `SpatRaster` objects. They must have the same extent and
resolution. However, if `x` is empty (has no cell values), its geometry
is ignored with a warning. Two empty SpatRasters with the same geometry
can also be combined (to get a summed number of layers). Also see
`add<-`

– Add a `SpatRaster` to a `SpatRasterDataset` or `SpatRasterCollection`

– Add `SpatVector` objects to a new or existing `SpatVectorCollection`

To append SpatVectors, use `rbind`.

## Usage

``` r
# S4 method for class 'SpatRaster'
c(x, ..., warn=TRUE)

# S4 method for class 'SpatRasterDataset'
c(x, ...)

# S4 method for class 'SpatRasterCollection'
c(x, ...)

# S4 method for class 'SpatVector'
c(x, ...)

# S4 method for class 'SpatVectorCollection'
c(x, ...)
```

## See also

`add<-`

## Arguments

- x:

  SpatRaster, SpatVector, SpatRasterDataset or SpatVectorCollection

- warn:

  logical. If `TRUE`, a warning is emitted if `x` is an empty SpatRaster

- ...:

  as for `x` (you can only combine raster with raster data and vector
  with vector data)

## Value

Same class as `x`

## Examples

``` r
r <- rast(nrows=5, ncols=9)
values(r) <- 1:ncell(r)
x <- c(r, r*2, r*3)
```
