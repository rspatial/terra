# Compare geometries

Evaluate whether two SpatRasters have the same extent, number of rows
and columns, projection, resolution, and origin (or a subset of these
comparisons).

Or evaluate whether two SpatVectors have the same geometries, or whether
a SpatVector has duplicated geometries.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
compareGeom(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE,
  rowcol=TRUE, res=FALSE, stopOnError=TRUE, messages=FALSE)

# S4 method for class 'SpatVector,SpatVector'
compareGeom(x, y, tolerance=0)

# S4 method for class 'SpatVector,missing'
compareGeom(x, y, tolerance=0)
```

## Arguments

- x:

  SpatRaster or SpatVector

- y:

  Same as `x`. If `x` is a SpatRaster, `y` can also be a list of
  SpatRasters. If `x` is a SpatVector, `y` can be missing

- ...:

  Additional SpatRasters

- lyrs:

  logical. If `TRUE`, the number of layers is compared

- crs:

  logical. If `TRUE`, coordinate reference systems are compared

- warncrs:

  logical. If `TRUE`, a warning is given if the crs is different
  (instead of an error)

- ext:

  logical. If `TRUE`, bounding boxes are compared

- rowcol:

  logical. If `TRUE`, number of rows and columns of the objects are
  compared

- res:

  logical. If `TRUE`, resolutions are compared (redundant when checking
  extent and rowcol)

- stopOnError:

  logical. If `TRUE`, code execution stops if raster do not match

- messages:

  logical. If `TRUE`, warning/error messages are printed even if
  `stopOnError=FALSE`

- tolerance:

  numeric

## Value

logical (SpatRaster) or matrix of logical (SpatVector)

## Examples

``` r
r1 <- rast()
r2 <- rast()
r3 <- rast()
compareGeom(r1, r2, r3)
#> [1] TRUE
nrow(r3) <- 10


if (FALSE) { # \dontrun{
compareGeom(r1, r3)
} # }
```
