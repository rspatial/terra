# Write raster data to a NetCDF file

Write a SpatRaster or SpatRasterDataset to a NetCDF file.

When using a SpatRasterDataset, the varname, longname, and unit should
be set in the object (see examples).

Always use the `".nc"` or `".cdf"` file extension to assure that the
file can be properly read again by GDAL

You can write multiple rasters (variables) that are two (x, y), three
(x, y, z or x, y, time) or four dimensional (x, y, z, time).

See [`depth`](https://rspatial.github.io/terra/reference/depth.md) and
[`time`](https://rspatial.github.io/terra/reference/time.md) for
specifying the axes of the thrid and/or fourth dimension(s).

## Usage

``` r
# S4 method for class 'SpatRaster'
writeCDF(x, filename, varname, longname="", unit="", split=FALSE, ...)

# S4 method for class 'SpatRasterDataset'
writeCDF(x, filename, overwrite=FALSE, timename="time", atts="", 
    gridmap="", prec="float", compression=NA, missval, tags=FALSE, ...)
```

## Arguments

- x:

  SpatRaster or SpatRasterDataset

- filename:

  character. Output filename

- varname:

  character. Name of the dataset

- longname:

  character. Long name of the dataset

- unit:

  character. Unit of the data

- split:

  logical. If `TRUE` each layer of `x` is treated as a sub-dataset

- atts:

  character. A vector of additional global attributes to write. The must
  be formatted like c("x=a value", "y=abc")

- gridmap:

  character. The crs is always written to the file in standard formats.
  With this argument you can also write the format commonly used in
  netcdf files. Something like
  `c("grid_mapping_name=lambert_azimuthal_equal_area", "longitude_of_projection_origin=10", "latitude_of_projection_origin=52", "false_easting=4321000", "false_northing=3210000")`

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- timename:

  character. The name of the "time" dimension

- prec:

  character. One of "double", "float", "integer", "short", "byte" or
  "char"

- compression:

  Can be set to an integer between 1 (least compression) and 9 (most
  compression)

- missval:

  numeric, the number used to indicate missing values

- tags:

  logical. If `TRUE` the value returned by
  [`metags`](https://rspatial.github.io/terra/reference/metags.md) are
  written to the file as attributes

- ...:

  additional arguments passed on to the SpatRasterDataset method, and
  from there possibly to
  [`ncvar_def`](https://rdrr.io/pkg/ncdf4/man/ncvar_def.html)

## Value

SpatRaster or SpatDataSet

## See also

see
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)
for writing other file formats

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
fname <- paste0(tempfile(), ".nc")
rr <- writeCDF(r, fname, overwrite=TRUE, varname="alt", 
      longname="elevation in m above sea level", unit="m")

a <- rast(ncols=5, nrows=5, nl=50)
values(a) <- 1:prod(dim(a))
time(a) <- as.Date("2020-12-31") + 1:nlyr(a)
aa <- writeCDF(a, fname, overwrite=TRUE, varname="power", 
      longname="my nice data", unit="U/Pa")

b <- sqrt(a)
s <- sds(a, b)
names(s) <- c("temp", "prec")
longnames(s) <- c("temperature (C)", "precipitation (mm)")
units(s) <- c("°C", "mm")
ss <- writeCDF(s, fname, overwrite=TRUE)

# four dimensional 
r1 <- rast(nrow=5, ncol=5, vals=1:100, nlyr=4)
depth(r1) <- c(0, 2, 0, 2)
time(r1) <- c(as.Date("2012-12-12") + c(1,1,2,2))
depthName(r1) <- "angle"

r2 <- rast(nrow=5, ncol=5, vals=1:150, nlyr=6)
depth(r2) <- c(10, 10, 20, 20, 30, 30)
time(r2) <- c(as.Date("2012-12-12") + c(1:2, 1:2, 1:2))
depthName(r2) <- "height"
depthUnit(r2) <- "cm"

s <- sds(r1, r2)
names(s) <- c("TH", "DBZH")
units(s) <- c("-", "Pa")
x <- writeCDF(s, filename = fname, overwrite=TRUE)
#> Warning: GDAL Message 1: dimension #1 (angle) is not a Time dimension.
#> Warning: GDAL Message 1: dimension #1 (time) is not a Time dimension.
#> Warning: GDAL Message 1: dimension #0 (height) is not a Time dimension.
x[1]
#> class       : SpatRaster 
#> size        : 5, 5, 4  (nrow, ncol, nlyr)
#> resolution  : 72, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source      : file209f56205edc.nc:TH 
#> varname     : TH 
#> names       : TH_angle=0_1, TH_angle=2_1, TH_angle=0_2, TH_angle=2_2 
#> unit        : - 
#> depth       : 0 to 2 (angle: 2 steps) 
#> time (days) : 2012-12-13 to 2012-12-14 (2 steps) 
time(x[1])
#> [1] "2012-12-13" "2012-12-13" "2012-12-14" "2012-12-14"
depth(x[1])
#> [1] 0 2 0 2

x[2]
#> class       : SpatRaster 
#> size        : 5, 5, 6  (nrow, ncol, nlyr)
#> resolution  : 72, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source      : file209f56205edc.nc:DBZH 
#> varname     : DBZH 
#> names       : DBZH_~=10_1, DBZH_~=10_2, DBZH_~=20_1, DBZH_~=20_2, DBZH_~=30_1, DBZH_~=30_2 
#> unit        : Pa 
#> depth       : 10 to 30 (height [cm]: 3 steps) 
#> time (days) : 2012-12-13 to 2012-12-14 (2 steps) 
time(x[2])
#> [1] "2012-12-13" "2012-12-14" "2012-12-13" "2012-12-14" "2012-12-13"
#> [6] "2012-12-14"
depth(x[2])
#> [1] 10 10 20 20 30 30

# for CRAN
file.remove(fname)
#> [1] TRUE
```
