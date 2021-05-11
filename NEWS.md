
# version 1.2-8

## new

- `as.lines` method for SpatRaster
- `as.polygons` method for SpatVector lines
- `autocor,numeric-method` has new methods `mean`, to compute the local mean, and `locmor`, for the local Moran's *I* 
- `sharedPaths` method for SpatVector (lines and polygons)
- `RGB2col` method to reduce a three-layer RGB SpatRaster to a single layer SpatRaster with a color-table (with <= 256 colors)
- `split` methods for SpatVector and SpatRaster

## enhancements

- `rast(Raster*)` now takes the crs from the Raster object, not from the file it may point to. Suggested by Floris Vanderhaeghe [#200](https://github.com/rspatial/terra/issues/200)
- `convhull` has a new argument `by=""` to make convex hulls for sub-sets of a SpatVector.


## bug fixes

- `crop` works again with `sf` objects. Reported by Sebastian Brinkmann [#201] (https://github.com/rspatial/terra/issues/201)
- `vect,sf-method` now also works for lines, and should be faster
- `vect,character` crashed R if a file had empty geometries. Reported by consumere [#202](https://github.com/rspatial/terra/issues/202)
- `extract(points, bilinear=TRUE, cells=TRUE)` now works. Reported by fab4app [#203](https://github.com/rspatial/terra/issues/203)


## name changes

To avoid name conflicts with the `spatstat` package

- `area,SpatRaster-method(x, sum=FALSE)` -> `cellSize(x)`
- `area,SpatRaster/SpatVector-method(x, sum=TRUE)` -> `expanse(x)`
- `convexhull` -> `convHull`
- `perimeter` -> `perim`
- `tiles` -> `makeTiles`
- `coords` -> `crds`


# version 1.2-5

## new

- `trim` has a new argument `value` that allows trimming rows and columns with other values than the default `NA`
- `rapp` has a new argument `clamp` that allows clamping start and end values to `1:nlyr(x)`, avoiding that all values are considered `NA`
- `spatSample,SpatRaster-method` has new arguments `as.points` and `values`. Getting values, cells and coordinates is no longer mutually exclusive. In response to [#191](https://github.com/rspatial/terra/issues/191). Requested by Agustin Lobo
- `area,SpatRaster-method` has a new argument `mask=FALSE`
- `classify` can now take a single number to request that many cuts
- `mosaic` and `merge` now warn and resample if rasters are not aligned
- `extract` has a new argument `exact` to get the fraction covered for each cell

## bug fixes

- `flip(x, direction="vertical")` no longer reverses the order of the layers
- `extract` did not work for horizontal or vertical lines as their extent was considered invalid. Reported by Monika Tomaszewska
- `autocor` did not handle NA values [#192](https://github.com/rspatial/terra/issues/192). Reported by Laurence Hawker
- `nearest` now works for angular coordinates
- The unit of `slope` in `terrain` was not correct (the tangent was returned instead of the slope) [#196](https://github.com/rspatial/terra/issues/196). Reported by Sven Alder
- `quantile` now works for rasters that have cells that are all `NA`. Reported by Jerry Nelson

## name changes

To avoid name conflicts with tidyverse 

with deprecation warning:

- separate -> segregate
- expand -> extend
- near -> nearby
- pack -> wrap 

without deprecation warning:

- transpose -> trans
- collapse -> tighten 
- fill -> fillHoles
- select -> sel


# version 1.1-17

## major changes 

- `c-SpatVector-method` now returns a list. `rbind` is used to append SpatVector objects
- overhaul of handling of factors. `rats` has been removed, and `levels` and `cats` have changed


# version 1.1-4

- No news recorded for this version and earlier versions
