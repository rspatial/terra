# SpatRaster class

A `SpatRaster` represents a rectangular part of the world that is
sub-divided into rectangular cells of equal area (in terms of the units
of the coordinate reference system). For each cell can have multiple
values ("layers").

An object of the `SpatRaster` class can point to one or more files on
disk that hold the cell values, and/or it can hold these values in
memory. These objects can be created with the
[`rast`](https://rspatial.github.io/terra/reference/rast.md) method.

A `SpatRasterDataset` is a collection of sub-datasets, where each is a
`SpatRaster` for the same area (extent) and coordinate reference system,
but possibly with a different resolution. Sub-datasets are often used to
capture variables (e.g. temperature and precipitation), or a fourth
dimension (e.g. height, depth or time) if the sub-datasets already have
three dimensions (multiple layers).

A `SpatRasterCollection` is a collection of SpatRasters with no
restriction in the extent or other geometric parameters.

## Examples

``` r
rast()
#> class       : SpatRaster 
#> size        : 180, 360, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
```
