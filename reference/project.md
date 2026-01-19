# Change the coordinate reference system

Change the coordinate reference system ("project") of a SpatVector,
SpatRaster or a matrix with coordinates.

## Usage

``` r
# S4 method for class 'SpatVector'
project(x, y, partial = FALSE)

# S4 method for class 'SpatRaster'
project(x, y, method, mask=FALSE, align_only=FALSE, res=NULL, 
  origin=NULL, threads=FALSE, filename="", ..., use_gdal=TRUE, by_util = FALSE)

# S4 method for class 'SpatExtent'
project(x, from, to)

# S4 method for class 'matrix'
project(x, from, to)
```

## Arguments

- x:

  SpatRaster, SpatVector, SpatExtent or matrix (with x and y columns)
  whose coordinates to project

- y:

  if `x` is a SpatRaster, the preferred approach is for `y` to be a
  SpatRaster as well, serving as a template for the geometry (extent and
  resolution) of the output SpatRaster. Alternatively, you can provide a
  coordinate reference system (CRS) description.

  You can use the following formats to define coordinate reference
  systems: WKT, PROJ.4 (e.g., `+proj=longlat +datum=WGS84`), or an EPSG
  code (e.g., `"epsg:4326"`). But note that the PROJ.4 notation has been
  deprecated, and you can only use it with the WGS84/NAD83 and NAD27
  datums. Other datums are silently ignored.

  If `x` is a SpatVector, you can provide a crs definition as discussed
  above, or any other object from which such a crs can be extracted with
  [`crs`](https://rspatial.github.io/terra/reference/crs.md)

- partial:

  logical. If `TRUE`, geometries that can only partially be represented
  in the output crs are included in the output

- method:

  character. Method used for estimating the new cell values of a
  SpatRaster. One of:

  `bilinear`: bilinear interpolation (3x3 cell window). This is used by
  default if the first layer of `x` is not categorical

  `mean`: This can be a good choice with continuous variables if the
  output cells overlap with multiple input cells.

  `near`: nearest neighbor. This is used by default if the first layer
  of `x` is categorical. This method is not a good choice for continuous
  values.

  `mode`: The modal value. This can be a good choice for categorical
  rasters, if the output cells overlap with multiple input cells.

  `cubic`: cubic interpolation (5x5 cell window).

  `cubicspline`: cubic B-spline interpolation. (5x5 cell window).

  `lanczos`: Lanczos windowed sinc resampling. (7x7 cell window).

  `sum`: the weighted sum of all non-NA contributing grid cells.

  `min, q1, median, q3, max`: the minimum, first quartile, median, third
  quartile, or maximum value.

  `rms`: the root-mean-square value of all non-NA contributing grid
  cells.

- mask:

  logical. If `TRUE`, mask out areas outside the input extent. For
  example, to avoid data wrapping around the date-line (see example with
  Robinson projection). To remove cells that are `NA` in `y` (if `y` is
  a SpatRaster) you can use the
  [`mask`](https://rspatial.github.io/terra/reference/mask.md)` method`
  after calling `project` (this function)

- align_only:

  logical. If `TRUE`, and `y` is a SpatRaster, the template is used for
  the spatial resolution and origin, but the extent is set such that all
  of the extent of `x` is included

- res:

  numeric. Can be used to set the resolution of the output raster if `y`
  is a CRS

- origin:

  numeric. Can be used to set the origin of the output raster if `y` is
  a CRS

- threads:

  logical. If `TRUE` multiple threads are used (faster for large files)

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

- use_gdal:

  logical. If `TRUE` the GDAL-warp algorithm is used. Otherwise, a
  slower internal algorithm is used that may be more accurate if there
  is much variation in the cell sizes of the output raster. Only the
  `near` and `bilinear` algorithms are available for the internal
  algorithm

- by_util:

  logical. If `TRUE` and `gdal=TRUE`, the GDAL warp utility is used

- from:

  character. Coordinate reference system of `x`

- to:

  character. Output coordinate reference system

## Value

SpatVector or SpatRaster

## See also

[`crs`](https://rspatial.github.io/terra/reference/crs.md),
[`resample`](https://rspatial.github.io/terra/reference/resample.md)

## Note

The PROJ.4 notation of coordinate reference systems has been partly
deprecated in the GDAL/PROJ library that is used by this function. You
can still use this notation, but \*only\* with the WGS84 datum. Other
datums are silently ignored.

Transforming (projecting) raster data is fundamentally different from
transforming vector data. Vector data can be transformed and
back-transformed without loss in precision and without changes in the
values. This is not the case with raster data. In each transformation
the values for the new cells are estimated in some fashion. Therefore,
if you need to match raster and vector data for analysis, you should
generally transform the vector data.

When using this method with a `SpatRaster`, the preferable approach is
to provide a template `SpatRaster` as argument `y`. The template is then
another raster dataset that you want your data to align with. If you do
not have a template to begin with, you can do `project(rast(x), crs)`
and then manipulate the output to get the template you want. For
example, where possible use whole numbers for the extent and resolution
so that you do not have to worry about small differences in the future.
You can use commands like `dim(z) = c(180, 360)` or `res(z) <- 100000`.

The output resolution should generally be similar to the input
resolution, but there is no "correct" resolution in raster
transformation. It is not obvious what this resolution is if you are
using lon/lat data that spans a large North-South extent.

## Examples

``` r
## SpatRaster
a <- rast(ncols=40, nrows=40, xmin=-110, xmax=-90, ymin=40, ymax=60, 
          crs="+proj=longlat +datum=WGS84")
values(a) <- 1:ncell(a)
newcrs="+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
b <- rast(ncols=94, nrows=124, xmin=-944881, xmax=935118, ymin=4664377, ymax=7144377, crs=newcrs)
w <- project(a, b)


## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
crs(v, proj=TRUE)
#> [1] "+proj=longlat +datum=WGS84 +no_defs"
cat(crs(v), "\n")
#> GEOGCRS["WGS 84",
#>     DATUM["World Geodetic System 1984",
#>         ELLIPSOID["WGS 84",6378137,298.257223563,
#>             LENGTHUNIT["metre",1]]],
#>     PRIMEM["Greenwich",0,
#>         ANGLEUNIT["degree",0.0174532925199433]],
#>     CS[ellipsoidal,2],
#>         AXIS["geodetic latitude (Lat)",north,
#>             ORDER[1],
#>             ANGLEUNIT["degree",0.0174532925199433]],
#>         AXIS["geodetic longitude (Lon)",east,
#>             ORDER[2],
#>             ANGLEUNIT["degree",0.0174532925199433]],
#>     ID["EPSG",4326]] 

project(v, "+proj=moll")
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 12, 6  (geometries, attributes)
#>  extent      : 437476.4, 497805.3, 5815524, 5892478  (xmin, xmax, ymin, ymax)
#>  coord. ref. : +proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04


project(v, "EPSG:2169")
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 12, 6  (geometries, attributes)
#>  extent      : 49540.31, 105922, 57009.53, 138631.1  (xmin, xmax, ymin, ymax)
#>  coord. ref. : LUREF / Luxembourg TM (EPSG:2169) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04
```
