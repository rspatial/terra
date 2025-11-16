# Create a SpatRaster

Methods to create a SpatRaster. These objects can be created from
scratch, from a filename, or from another object.

A SpatRaster represents a spatially referenced surface divided into
three dimensional cells (rows, columns, and layers).

When a SpatRaster is created from one or more files, it does not load
the cell (pixel) values into memory (RAM). It only reads the parameters
that describe the geometry of the SpatRaster, such as the number of rows
and columns and the coordinate reference system. The actual values will
be read when needed.

Note that there are operating system level limitations to the number of
files that can be opened simultaneously. Using a SpatRaster of very many
files (e.g. 10,000) may cause R to crash when you use it in a
computation. In situations like that you may need to split up the task
or combine data into fewer (multi-layer) files. Also note that the GTiff
format used for temporary files cannot store more than 65,535 layers in
a single file.

## Usage

``` r
# S4 method for class 'character'
rast(x, subds=0, lyrs=NULL, drivers=NULL, opts=NULL, win=NULL, 
    snap="near", vsi=FALSE, raw=FALSE, noflip=FALSE, 
    guessCRS=TRUE, domains="", md=FALSE, dims=NULL)

# S4 method for class 'missing'
rast(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, ymin=-90,
    ymax=90, crs, extent, resolution, vals, names, time, units)

# S4 method for class 'SpatRaster'
rast(x, nlyrs=nlyr(x), names, vals, keeptime=TRUE, 
    keepunits=FALSE, props=FALSE, tags=FALSE) 

# S4 method for class 'matrix'
rast(x, type="", crs="", digits=6, extent=NULL)

# S4 method for class 'data.frame'
rast(x, type="xyz", crs="", digits=6, extent=NULL)

# S4 method for class 'array'
rast(x, crs="", extent=NULL)

# S4 method for class 'list'
rast(x, warn=TRUE)

# S4 method for class 'SpatRasterDataset'
rast(x)

# S4 method for class 'SpatVector'
rast(x, type="", ...)
                    
# S4 method for class 'SpatExtent'
rast(x, ...)
```

## Arguments

- x:

  filename (character), missing, SpatRaster, SpatRasterDataset,
  SpatExtent, SpatVector, matrix, array, list of SpatRasters. For other
  types it will be attempted to create a SpatRaster via (\`as(x,
  "SpatRaster")\`

- subds:

  positive integer or character to select a sub-dataset. If zero or "",
  all sub-datasets are returned (if possible)

- lyrs:

  positive integer or character to select a subset of layers (a.k.a.
  "bands"). If `x` has multiple filenames, the same layer numbers are
  selected from each of the files, unless numbers larger than the number
  of layers of the first data source are included

- drivers:

  character. GDAL drivers to consider

- opts:

  character. GDAL dataset open options

- win:

  SpatExtent to set a
  [`window`](https://rspatial.github.io/terra/reference/window.md) (area
  of interest)

- snap:

  character. One of "near", "in", or "out", to indicate how the extent
  of [`window`](https://rspatial.github.io/terra/reference/window.md)
  should be "snapped" to `x`

- vsi:

  logical. If `TRUE`, "\vsicurl\\ is prepended to filenames that start
  with "http". There are many [VSI configuration
  options](https://gdal.org/en/stable/user/virtual_file_systems.html)
  that can be set with
  [`setGDALconfig`](https://rspatial.github.io/terra/reference/gdal.md)

- raw:

  logical. If `TRUE`, scale and offset values are ignored. See
  [`scoff`](https://rspatial.github.io/terra/reference/scoff.md) to get
  these parameters

- noflip:

  logical. If `TRUE`, a raster (e.g. JPEG image) that is not
  georeferenced and that GDAL assigns a flipped extent to
  (`ymax < ymin`), is not considered flipped. This avoids the need to
  [`flip`](https://rspatial.github.io/terra/reference/flip.md) the
  raster vertically

- guessCRS:

  logical. If `TRUE` and the the file does not specify a CRS but has an
  extent that is within longitude/latitude bounds, the
  longitude/latitude crs is assigned to the SpatRaster

- domains:

  character. Metadata domains to read (see
  [`metags`](https://rspatial.github.io/terra/reference/metags.md) to
  retrieve their values if there are any. "" is the default domain

- md:

  logical. If `TRUE`, the multi-dimensional GDAL interface is used under
  the hood for file reading. This interface can only be used for a few
  file formats (netCDF/HDF5) and can sometimes (not always) provide
  notably faster reading of data with many (time) steps in the third or
  higher dimension. Support for this is new and experimental (June 2025)

- dims:

  numeric. Specify the order of the dimensions to read atypical files.
  See
  [`ar_info`](https://rspatial.github.io/terra/reference/ar_info.md).
  Only relevant if `md=TRUE`. Not used yet

- nrows:

  positive integer. Number of rows

- ncols:

  positive integer. Number of columns

- nlyrs:

  positive integer. Number of layers

- xmin:

  minimum x coordinate (left border)

- xmax:

  maximum x coordinate (right border)

- ymin:

  minimum y coordinate (bottom border)

- ymax:

  maximum y coordinate (top border)

- crs:

  character. Description of the Coordinate Reference System (map
  projection) in `PROJ.4`, `WKT` or `authority:code` notation. See
  [`crs`](https://rspatial.github.io/terra/reference/crs.md). If this
  argument is missing, and the x coordinates are within -360 .. 360 and
  the y coordinates are within -90 .. 90, longitude/latitude is assigned

- keeptime:

  logical. If `FALSE` the time stamps are discarded

- keepunits:

  logical. If `FALSE` the layer units are discarded

- props:

  logical. If `TRUE` the properties (categories and color-table) are
  kept

- tags:

  logical. If `TRUE` the user specified metadata tags are kept (see
  [`metags`](https://rspatial.github.io/terra/reference/metags.md)).

- extent:

  object of class SpatExtent. If present, the arguments xmin, xmax, ymin
  and ymax are ignored

- resolution:

  numeric vector of length 1 or 2 to set the spatial resolution (see
  [`res`](https://rspatial.github.io/terra/reference/dimensions.md)). If
  this argument is used, arguments `ncols` and `nrows` are ignored

- vals:

  numeric. An optional vector with cell values (if fewer values are
  provided, these are recycled to reach the number of cells)

- names:

  character. An optional vector with layer names (must match the number
  of layers)

- time:

  time or date stamps for each layer

- units:

  character. units for each layer

- type:

  character. If the value is `"xyz"`, `x` must be a SpatVector with
  point geometry, or a matrix or data.frame with at least two columns,
  the first with `x` (or longitude) and the second with `y` (or
  latitude) coordinates that represent the centers of raster cells. The
  additional columns are the values associated with the raster cells. If
  the value is `"xylz"`, `x` must have four columns with the third
  representing the layer and the fourth the cell values. If the value is
  `""`, the resulting SpatRaster will have the same number of rows and
  columns as `x`.

- digits:

  integer to set the precision for detecting whether points are on a
  regular grid (a low number of digits is a low precision). Only used
  when `type="xyz"`

- warn:

  logical. If `TRUE`, a warnings about empty rasters may be emitted

- ...:

  additional arguments passed on to the `rast,missing-method`

## Value

SpatRaster

## Details

Files are read with the GDAL library. GDAL guesses the file format from
the name, and/or tries reading it with different "drivers" (see
[`gdal`](https://rspatial.github.io/terra/reference/gdal.md)) until it
succeeds. In very few cases this may cause a file to be opened with the
wrong driver, and some information may be lost. For example, when a
netCDF file is opened with the HDF5 driver. You can avoid that by using
argument `rast("filename.ncdf", drivers="NETCDF")`

These classes hold a C++ pointer to the data "reference class" and that
creates some limitations. They cannot be recovered from a saved R
session either or directly passed to nodes on a computer cluster.
Generally, you should use
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)
to save SpatRaster objects to disk (and pass a filename or cell values
of cluster nodes). Also see
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md).

## See also

[`sds`](https://rspatial.github.io/terra/reference/sds.md) to create a
SpatRasterDataset (SpatRasters with the same geometry representing
different variables or higher dimension),
[`sprc`](https://rspatial.github.io/terra/reference/sprc.md) to create a
SpatRasterCollection (to combine SpatRasters with different geometries),
and [`vect`](https://rspatial.github.io/terra/reference/vect.md) for
vector (points, lines, polygons) data

## Examples

``` r
# Create a SpatRaster from scratch
x <- rast(nrows=108, ncols=21, xmin=0, xmax=10)

# Create a SpatRaster from a file
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

# A file with multiple layers. This one is special as the layers are RGB color channels 
s <- rast(system.file("ex/logo.tif", package="terra"))

# remove the color channels
#plot(s)
#RGB(s) <- NULL
#plot(s)

# Create a skeleton with no associated cell values
rast(s)
#> class       : SpatRaster 
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 

# from a matrix 
m <- matrix(1:25, nrow=5, ncol=5)
rm <- rast(m)

# from a "xyz" data.frame
d <- as.data.frame(rm, xy=TRUE)
head(d)
#>     x   y lyr.1
#> 1 0.5 4.5     1
#> 2 1.5 4.5     6
#> 3 2.5 4.5    11
#> 4 3.5 4.5    16
#> 5 4.5 4.5    21
#> 6 0.5 3.5     2
rast(d, type="xyz")
#> class       : SpatRaster 
#> size        : 5, 5, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 5, 0, 5  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> name        : lyr.1 
#> min value   :     1 
#> max value   :    25 
```
