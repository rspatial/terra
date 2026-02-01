# Make tiles or get their extents

Divide a SpatRaster into "tiles". The cells of another SpatRaster
(normally with a much lower resolution) or a SpatVector with polygon
geometry can be used to define the tiles. You can also provide one or
two numbers to indicate the number of rows and columns per tile.

`getTileExtents` returns the extents of the (virtual) tiles, while
`makeTiles` creates files for the tiles and returns their filenames.

## Usage

``` r
# S4 method for class 'SpatRaster'
makeTiles(x, y, filename="tile_.tif", extend=FALSE,
    na.rm=FALSE, buffer=0, value="files", overwrite=FALSE, ...)


# S4 method for class 'SpatRaster'
getTileExtents(x, y, extend=FALSE, buffer=0)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatRaster or SpatVector defining the zones; or a positive integer
  specifying the number of rows and columns for each zone (or 2 numbers
  to differentiate the number of rows and columns)

- filename:

  character. Output filename template. Filenames will be altered by
  adding the tile number for each tile

- extend:

  logical. If `TRUE`, the extent of `y` is expanded to assure that it
  covers all of `x`

- na.rm:

  logical. If `TRUE`, tiles with only missing values are ignored

- buffer:

  integer. The number of additional rows and columns added to each tile.
  Can be a single number, or two numbers to specify a separate number of
  rows and columns. This allows for creating overlapping tiles that can
  be used for computing spatial context dependent values with e.g.
  [`focal`](https://rspatial.github.io/terra/reference/focal.md). The
  expansion is only inside `x`, no rows or columns outside of `x` are
  added

- value:

  character. The type of return value desired. Either "files" (for the
  filenames), "raster" (for a SpatRaster), or "collection" (for a
  SpatRasterCollection)

- overwrite:

  logical. If `TRUE`, existing tiles are overwritten; otherwise they are
  skipped (without error or warning)

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

`makeTiles` returns a character (filenames), SpatRaster or
SpatRasterCollection value. `getTileExtents` returns a matrix with
extents

## See also

[`vrt`](https://rspatial.github.io/terra/reference/vrt.md) to create a
SpatRaster from tiles;
[`crop`](https://rspatial.github.io/terra/reference/crop.md) for
sub-setting arbitrary parts of a SpatRaster;
[`divide`](https://rspatial.github.io/terra/reference/divide.md) to
divide a SpatRaster into parts.

## Examples

``` r
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)

getTileExtents(r, x)
#>      xmin xmax ymin ymax
#> [1,] -180    0    0   90
#> [2,]    0  180    0   90
#> [3,] -180    0  -90    0
#> [4,]    0  180  -90    0
getTileExtents(r, x, buffer=3)
#>        xmin  xmax  ymin ymax
#> [1,] -190.8  10.8  -5.4 95.4
#> [2,]  -10.8 190.8  -5.4 95.4
#> [3,] -190.8  10.8 -95.4  5.4
#> [4,]  -10.8 190.8 -95.4  5.4


filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff
#> [1] "/tmp/Rtmp20FJ70/file236a1704533c_1.tif"
#> [2] "/tmp/Rtmp20FJ70/file236a1704533c_2.tif"
#> [3] "/tmp/Rtmp20FJ70/file236a1704533c_3.tif"
#> [4] "/tmp/Rtmp20FJ70/file236a1704533c_4.tif"

vrt(ff)
#> class       : SpatRaster 
#> size        : 100, 100, 1  (nrow, ncol, nlyr)
#> resolution  : 3.6, 1.8  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : spat_236a33c025be_9066_DvggOkjkTOQeAAI.vrt 
#> name        : spat_236a33c025be_9066_DvggOkjkTOQeAAI 
#> min value   :                                      1 
#> max value   :                                  10000 
```
