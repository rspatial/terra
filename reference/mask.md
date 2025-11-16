# Mask values in a SpatRaster or SpatVector

If `x` is a `SpatRaster`: Create a new SpatRaster that has the same
values as SpatRaster `x`, except for the cells that are `NA` (or other
`maskvalue`) in another SpatRaster (the 'mask'), or the cells that are
not covered by a SpatVector or SpatExtent. These cells become `NA` (or
another `updatevalue`).

If `x` is a SpatVector or SpatExtent: Select geometries of `x` that
intersect, or not intersect, with the geometries of `y`.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
mask(x, mask, inverse=FALSE, maskvalues=NA, 
   updatevalue=NA, filename="", ...)
   
# S4 method for class 'SpatRaster,SpatVector'
mask(x, mask, inverse=FALSE, updatevalue=NA,
  touches=TRUE, filename="", ...)

# S4 method for class 'SpatRaster,SpatExtent'
mask(x, mask, inverse=FALSE, updatevalue=NA,
  touches=TRUE, filename="", ...)

# S4 method for class 'SpatVector,SpatVector'
mask(x, mask, inverse=FALSE)

# S4 method for class 'SpatVector,SpatExtent'
mask(x, mask, inverse=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- mask:

  SpatRaster or SpatVector

- inverse:

  logical. If `TRUE`, areas on mask that are \_not\_ the `maskvalue` are
  masked

- maskvalues:

  numeric. The value(s) in `mask` that indicate which cells of `x`
  should be masked (change their value to `updatevalue` (default =
  `NA`))

- updatevalue:

  numeric. The value that masked cells should become (if they are not
  `NA`)

- touches:

  logical. If `TRUE`, all cells touched by lines or polygons will be
  masked, not just those on the line render path, or whose center point
  is within the polygon

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`subst`](https://rspatial.github.io/terra/reference/subst.md),
[`crop`](https://rspatial.github.io/terra/reference/crop.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
msk <- ifel(r < 400, NA, 1)

m <- mask(r, msk)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)[1,]

mv1 <- mask(r, v)
mv2 <- crop(r, v, mask=TRUE)
```
