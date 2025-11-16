# Merge SpatRasters, or merge a SpatVector with a data.frame

Merge multiple SpatRasters to create a new SpatRaster with a larger
spatial extent. The SpatRasters should all have the same coordinate
reference system. They should normally also have the same spatial origin
and resolution, but automatic resampling can be done depending on the
algorithm used (see argument `algo`). In areas where the SpatRasters
overlap, the values of the SpatRaster that is first in the sequence of
arguments (or in the SpatRasterCollection) will be retained (unless
`first=FALSE`).

There is also a method for merging SpatVector with a data.frame; that
is, to join the data.frame to the attribute table of the SpatVector.

See [`classify`](https://rspatial.github.io/terra/reference/classify.md)
to merge a SpatRaster with a `data.frame`.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
merge(x, y, ..., first=TRUE, na.rm=TRUE, algo=1, method=NULL, 
      filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRasterCollection,missing'
merge(x, first=TRUE, na.rm=TRUE, algo=1, method=NULL, filename="", ...)

# S4 method for class 'SpatVector,data.frame'
merge(x, y, ...)
```

## Arguments

- x:

  SpatRaster, SpatRasterCollection, or SpatVector

- y:

  missing if `x` is a SpatRasterCollection. SpatRaster if `x` is a
  SpatRaster. data.frame if `x` is a SpatVector

- ...:

  if `x` is a SpatRaster: additional objects of the same class as `x`.
  If `x` is a SpatRasterCollection: options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md).
  If `x` is a SpatVector, the same arguments as in
  [`merge`](https://rdrr.io/r/base/merge.html)

- first:

  logical. If `TRUE`, in areas where rasters overlap, the first value is
  used. Otherwise the last value is used

- na.rm:

  logical. If `TRUE` missing values are are ignored. This is only used
  for algo 1; the other two always ignore missing values

- algo:

  integer. You can use 1, 2 or 3 to pick a merge algorithm. algo 1 is
  generally faster than algo 2, but it may have poorer file compression.
  Algo 1 will resample input rasters (and that may slow it down), but
  algo 2 does not do that. You can increase the tolerance option to
  effectively get nearest neighbor resampling with, for example,
  `wopt=list(tolerance=.2)` allows misalignment of .2 times the
  resolution of the first input raster and effectively use nearest
  neighbor resampling. Algo 3 creates a virtual raster (see
  [`vrt`](https://rspatial.github.io/terra/reference/vrt.md)). This is
  very quick and can be a good approach if the merge raster is used as
  input to a next step in the analysis. It allows any amount of
  misalignment (and does not respond to the tolerance option). Otherwise
  its speed is similar to that of algo 2

- method:

  character. The interpolation method to be used if resampling is
  necessary (see argument `algo`). One of "nearest", "bilinear",
  "cubic", "cubicspline", "lanczos", "average", "mode" as in
  [`resample`](https://rspatial.github.io/terra/reference/resample.md).
  If `NULL`, "nearest" is used for categorical rasters and "bilinear"
  for other rasters

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster or SpatVector

## See also

Combining tiles with
[`vrt`](https://rspatial.github.io/terra/reference/vrt.md) may be more
efficient than using `merge`. See
[`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md) for
averaging overlapping regions.

See [`classify`](https://rspatial.github.io/terra/reference/classify.md)
to merge a `SpatRaster` and a `data.frame` and
[`union`](https://rspatial.github.io/terra/reference/union.md) to
combine SpatExtent objects.

## Examples

``` r
x <- rast(xmin=-110, xmax=-80, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-85, xmax=-55, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-60, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- merge(x, y, z)
m2 <- merge(z, y, x)
m3 <- merge(y, x, z)
# panel(c(m1, m2, m3))

# if you have many SpatRasters, it may be convenient
# to make a SpatRasterCollection
# s <- sprc(list(x, y, z))
s <- sprc(x, y, z)

sm1 <- merge(s, algo=1, first=FALSE)
sm2 <- merge(s, algo=2, first=FALSE)
#sm3 <- merge(s, algo=3, first=FALSE)

## SpatVector with data.frame
f <- system.file("ex/lux.shp", package="terra")
p <- vect(f)
dfr <- data.frame(District=p$NAME_1, Canton=p$NAME_2, Value=round(runif(length(p), 100, 1000)))
dfr <- dfr[1:5, ]
pm <- merge(p, dfr, all.x=TRUE, by.x=c('NAME_1', 'NAME_2'), by.y=c('District', 'Canton'))
pm
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 12, 7  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.44781, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1   NAME_2  ID_1  ID_2  AREA       POP Value
#>  type        :    <chr>    <chr> <num> <num> <num>     <num> <num>
#>  values      : Diekirch Clervaux     1     1   312 1.808e+04   796
#>                Diekirch Diekirch     1     2   218 3.254e+04   298
#>                Diekirch  Redange     1     3   259 1.866e+04   744
values(pm)
#>          NAME_1           NAME_2 ID_1 ID_2 AREA    POP Value
#> 1      Diekirch         Clervaux    1    1  312  18081   796
#> 2      Diekirch         Diekirch    1    2  218  32543   298
#> 3      Diekirch          Redange    1    3  259  18664   744
#> 4      Diekirch          Vianden    1    4   76   5163   698
#> 5      Diekirch            Wiltz    1    5  263  16735   736
#> 6  Grevenmacher       Echternach    2    6  188  18899    NA
#> 7  Grevenmacher           Remich    2    7  129  22366    NA
#> 8  Grevenmacher     Grevenmacher    2   12  210  29828    NA
#> 9    Luxembourg         Capellen    3    8  185  48187    NA
#> 10   Luxembourg Esch-sur-Alzette    3    9  251 176820    NA
#> 11   Luxembourg       Luxembourg    3   10  237 182607    NA
#> 12   Luxembourg           Mersch    3   11  233  32112    NA
```
