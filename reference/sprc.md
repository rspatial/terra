# Create a SpatRasterCollection

Methods to create a SpatRasterCollection. This is an object to hold a
collection (list) of SpatRasters. There are no restrictions on the
similarity of the SpatRaster geometry.

They can be used to combine several SpatRasters to be used with
[`merge`](https://rspatial.github.io/terra/reference/merge.md) or
[`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md)

You can create a SpatRasterCollection from a file with subdatasets.

## Usage

``` r
# S4 method for class 'character'
sprc(x, ids=0, opts=NULL, raw=FALSE, noflip=FALSE, guessCRS=TRUE, domains="") 

# S4 method for class 'SpatRaster'
sprc(x, ...) 

# S4 method for class 'list'
sprc(x) 

# S4 method for class 'missing'
sprc(x)
```

## Arguments

- x:

  SpatRaster, list with SpatRasters, missing, or filename

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

  logical. If `TRUE` and the the file does not specify a CRS but has an
  extent that is within longitude/latitude bounds, the
  longitude/latitude crs is assigned to the SpatRaster

- domains:

  character. Metadata domains to read (see
  [`metags`](https://rspatial.github.io/terra/reference/metags.md) to
  retrieve their values if there are any. "" is the default domain

- ...:

  additional SpatRasters

## Value

SpatRasterCollection

## See also

[`sds`](https://rspatial.github.io/terra/reference/sds.md)

## Examples

``` r
x <- rast(xmin=-110, xmax=-50, ymin=40, ymax=70, ncols=60, nrows=30)
y <- rast(xmin=-80, xmax=-20, ymax=60, ymin=30)
res(y) <- res(x)
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)

z <- sprc(x, y)
z
#> class       : SpatRasterCollection 
#> length      : 2 
#> nrow        : 30, 30 
#> ncol        : 60, 60 
#> nlyr        :  1,  1 
#> extent      : -110, -20, 30, 70  (xmin, xmax, ymin, ymax)
#> crs (first) : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
```
