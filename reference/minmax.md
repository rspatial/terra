# Get or compute the minimum and maximum cell values

The minimum and maximum value of a SpatRaster are returned or computed
(from a file on disk if necessary) and stored in the object.

## Usage

``` r
# S4 method for class 'SpatRaster'
minmax(x, compute=FALSE)
# S4 method for class 'SpatRaster'
hasMinMax(x)
# S4 method for class 'SpatRaster'
setMinMax(x, force=FALSE)
```

## Arguments

- x:

  SpatRaster

- compute:

  logical. If `TRUE` min and max values are computed if they are not
  available

- force:

  logical. If `TRUE` min and max values are recomputed even if already
  available

## Value

minmax: numeric matrix of minimum and maximum cell values by layer

hasMinMax: logical indicating whether the min and max values are
available.

setMinMax: nothing. Used for the side-effect of computing the minimum
and maximum values of a SpatRaster

## See also

[`where.min`](https://rspatial.github.io/terra/reference/where.md),
[`where.max`](https://rspatial.github.io/terra/reference/where.md)

## Examples

``` r
r <- rast(system.file("ex/elev.tif", package="terra"))
minmax(r)
#>     elevation
#> min       141
#> max       547
```
