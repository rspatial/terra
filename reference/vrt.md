# Virtual Raster Dataset

Create a Virtual Raster Dataset (VRT) from a collection of file-based
raster datasets (tiles). See
[gdalbuildvrt](https://gdal.org/en/latest/programs/gdalbuildvrt.html)
for details.

## Usage

``` r
# S4 method for class 'character'
vrt(x, filename="", options=NULL, overwrite=FALSE, set_names=FALSE, return_filename=FALSE)

# S4 method for class 'SpatRasterCollection'
vrt(x, filename="", options=NULL, overwrite=FALSE, return_filename=FALSE)
```

## Arguments

- x:

  SpatRasterCollection or character vector with filenames of raster
  "tiles". That is, files that have data for, typically non-overlapping,
  sub-regions of an raster. See
  [`makeTiles`](https://rspatial.github.io/terra/reference/makeTiles.md)

- filename:

  character. output VRT filename

- options:

  character. All arguments as separate vector elements. Options as for
  [gdalbuildvrt](https://gdal.org/en/latest/programs/gdalbuildvrt.html)

- overwrite:

  logical. Should `filename` be overwritten if it exists?

- set_names:

  logical. Add the layer names of the first tile to the vrt? If
  `options` includes `"-separate"` the name of each source file is
  added, and each input goes into a separate band in the VRT dataset

- return_filename:

  logical. If `TRUE` the filename is returned, otherwise a SpatRaster is
  returned

## Value

SpatRaster

## Note

A VRT can reference very many datasets. These are not all opened at the
same time. The default is to open not more than 100 files. To increase
performance, this maximum limit can be increased by setting the
GDAL_MAX_DATASET_POOL_SIZE configuration option to a bigger value with
[`setGDALconfig`](https://rspatial.github.io/terra/reference/gdal.md).
Note that a typical user process on Linux is limited to 1024
simultaneously opened files.

## See also

[`makeTiles`](https://rspatial.github.io/terra/reference/makeTiles.md)
to create tiles;
[`makeVRT`](https://rspatial.github.io/terra/reference/makeVRT.md) to
create a .vrt file for a binary raster file that does not have a header
file.
[`vrt_tiles`](https://rspatial.github.io/terra/reference/vrt_tiles.md)
to get the filenames of the tiles in a VRT.

## Examples

``` r
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)
filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff
#> [1] "/tmp/RtmpwQn5lJ/file2039521817a7_1.tif"
#> [2] "/tmp/RtmpwQn5lJ/file2039521817a7_2.tif"
#> [3] "/tmp/RtmpwQn5lJ/file2039521817a7_3.tif"
#> [4] "/tmp/RtmpwQn5lJ/file2039521817a7_4.tif"

#vrtfile <- paste0(tempfile(), ".vrt")
#v <- vrt(ff, vrtfile)


## output in lower resolution
#vrtfile <- paste0(tempfile(), ".vrt")
#v <- vrt(ff, vrtfile, options = c("-tr", 5, 5))
#head(readLines(vrtfile))
#v
```
