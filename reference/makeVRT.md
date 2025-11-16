# Make a VRT header file

Create a VRT header file for a "flat binary" raster file that needs a
header file to be able to read it, but does not have it.

## Usage

``` r
makeVRT(filename, nrow, ncol, nlyr=1, extent, xmin, ymin, xres, yres=xres, xycenter=TRUE,
   crs="+proj=longlat", lyrnms="", datatype, NAflag=NA, bandorder="BIL", byteorder="LSB",
   toptobottom=TRUE, offset=0, scale=1)
```

## Arguments

- filename:

  character. raster filename (without the ".vrt" extension)

- nrow:

  positive integer, the number of rows

- ncol:

  positive integer, the number of columns

- nlyr:

  positive integer, the number of layers

- extent:

  SpatExtent or missing

- xmin:

  numeric. minimum x coordinate (only used if `extent` is missing)

- ymin:

  numeric. minimum y coordinate (only used if `extent` is missing)

- xres:

  positive number. x resolution

- yres:

  positive number. y resolution)

- xycenter:

  logical. If `TRUE`, `xmin` and `xmax` represent the coordinates of the
  center of the extreme cell, in stead of the coordinates of the outside
  corner. Only used of `extent` is missing

- crs:

  character. Coordinate reference system description

- lyrnms:

  character. Layer names

- datatype:

  character. One of "INT2S", "INT4S", "INT1U", "INT2U", "INT4U",
  "FLT4S", "FLT8S". If missing, this is guessed from the file size
  (INT1U for 1 byte per value, INT2S for 2 bytes and FLT4S for 4 bytes
  per value). This may be wrong because, for example, 2 bytes per value
  may in fact be INT2U (with the U for unsigned) values

- NAflag:

  numeric. The value used as the "NA flag"

- bandorder:

  character. One of "BIL", "BIP", or "BSQ". That is Band Interleaved by
  Line, or by Pixel, or Band SeQuential

- byteorder:

  character. One of "LSB", "MSB". "MSB" is common for files generated on
  Linux systems, whereas "LSB" is common for files generated on windows

- toptobottom:

  logical. If `FALSE`, the values are read bottom to top

- offset:

  numeric. offset to be applied

- scale:

  numeric. scale to be applied

## Value

character (.VRT filename)

## See also

[`vrt`](https://rspatial.github.io/terra/reference/vrt.md) to create a
vrt for a collection of raster tiles
