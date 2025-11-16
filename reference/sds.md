# Create a SpatRasterDataset

Methods to create a SpatRasterDataset. This is an object to hold
"sub-datasets", each represented by a SpatRaster that may have multiple
layers. All sub-datasets must have the same raster geometry (extent and
resolution). You can use a SpatRasterCollection (see
[`sprc`](https://rspatial.github.io/terra/reference/sprc.md)) to combine
SpatRasters with different geometries.

See [`describe`](https://rspatial.github.io/terra/reference/describe.md)
for getting information about the sub-datasets present in a file.

## Usage

``` r
# S4 method for class 'missing'
sds(x) 

# S4 method for class 'character'
sds(x, ids=0, opts=NULL, raw=FALSE, noflip=FALSE, guessCRS=TRUE, domains="")

# S4 method for class 'SpatRaster'
sds(x, ...) 

# S4 method for class 'list'
sds(x) 

# S4 method for class 'array'
sds(x, crs="", extent=NULL)
```

## Arguments

- x:

  character (filename), or SpatRaster, or list of SpatRasters, or
  missing. If multiple filenames are provided, it is attempted to make
  SpatRasters from these, and combine them into a SpatRasterDataset

- ids:

  optional. vector of integer subdataset ids. Ignored if the first value
  is not a positive integer

- opts:

  character. GDAL dataset open options

- raw:

  logical. If `TRUE`, scale and offset values are ignored

- noflip:

  logical. If `TRUE`, a raster (e.g. JPEG image) that is not
  georeferenced and that GDAL assigns a flipped extent to
  (`ymax < ymin`), is not considered flipped. This avoids the need to
  [`flip`](https://rspatial.github.io/terra/reference/flip.md) the
  raster vertically

- guessCRS:

  logical. If `TRUE` and the file does not specify a CRS but has an
  extent that is within longitude/latitude bounds, the
  longitude/latitude crs is assigned to the SpatRaster

- domains:

  character. Metadata domains to read (see
  [`metags`](https://rspatial.github.io/terra/reference/metags.md) to
  retrieve their values if there are any). "" is the default domain

- crs:

  character. Description of the Coordinate Reference System (map
  projection) in `PROJ.4`, `WKT` or `authority:code` notation. If this
  argument is missing, and the x coordinates are within -360 .. 360 and
  the y coordinates are within -90 .. 90, longitude/latitude is assigned

- extent:

  [`SpatExtent`](https://rspatial.github.io/terra/reference/SpatExtent-class.md)

- ...:

  additional `SpatRaster` objects

## Value

SpatRasterDataset

## See also

[`sprc`](https://rspatial.github.io/terra/reference/sprc.md),
[`describe`](https://rspatial.github.io/terra/reference/describe.md)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- sds(s, s/2)
names(x) <- c("first", "second")
x
#> class       : SpatRasterDataset 
#> subdatasets : 2 
#> dimensions  : 77, 101 (nrow, ncol)
#> nlyr        : 3, 3 
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : logo.tif, memory 
#> names       : first, second 
length(x)
#> [1] 2

# extract the second SpatRaster
x[2]
#> class       : SpatRaster 
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> varname     : logo 
#> names       :   red, green,  blue 
#> min values  :   0.0,   0.0,   0.0 
#> max values  : 127.5, 127.5, 127.5 

a <- array(1:9, c(3,3,3,3))
sds(a)
#> class       : SpatRasterDataset 
#> subdatasets : 3 
#> dimensions  : 3, 3 (nrow, ncol)
#> nlyr        : 3, 3, 3 
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 3, 0, 3  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory 
```
