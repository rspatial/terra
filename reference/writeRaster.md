# Write raster data to a file

Write a SpatRaster to a file.

## Usage

``` r
# S4 method for class 'SpatRaster,character'
writeRaster(x, filename, overwrite=FALSE, ...)
```

## Arguments

- x:

  SpatRaster

- filename:

  character. Output filename. Can be a single filename, or as many
  filenames as `nlyr(x)` to write a file for each layer

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  additional arguments for writing files. See Details

## Value

SpatRaster. This function is used for the side-effect of writing values
to a file.

## See also

see [`writeCDF`](https://rspatial.github.io/terra/reference/writeCDF.md)
for writing NetCDF files.

## Details

In writeRaster, and in other methods that generate SpatRasters, options
for writing raster files to disk can be provided as additional arguments
or, in a few cases, as the `wopt` argument (a named list) if the
additional arguments are already used for a different purpose. See
[`terraOptions`](https://rspatial.github.io/terra/reference/terraOptions.md)
to get or set default values. The following options are available:

|            |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **name**   | **description**                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `datatype` | values accepted are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". With GDAL \>= 3.5 you can also use "INT8U" and "INT8S". And with GDAL \>= 3.7 you can use also use "INT1S". See [`gdal`](https://rspatial.github.io/terra/reference/gdal.md) to discover the GDAL version you are using. The first three letters indicate whether the datatype is an integer (whole numbers) of a real number ("float", decimal numbers), the fourth character indicates the number of bytes used for each number. Higher values allow for storing larger numbers and/or more precision; but create larger files. The "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (zero and positive values only). |
| `filetype` | file format expresses as [GDAL driver names](https://gdal.org/en/latest/drivers/raster/index.html). If this argument is not supplied, the driver is derived from the filename. You can use `gdal(drivers=TRUE)` to see what drivers are available in your installation                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `gdal`     | GDAL driver specific datasource creation options. See the GDAL documentation. For example, with the [GeoTiff file format](https://gdal.org/en/latest/drivers/raster/gtiff.html) you can use `gdal=c("COMPRESS=DEFLATE", "TFW=YES")`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `tempdir`  | the path where temporary files are to be written to.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `progress` | positive integer. If the number of chunks is larger, a progress bar is shown.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `memfrac`  | numeric between 0 and 0.9 (higher values give a warning). The fraction of available RAM that terra is allowed to use.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `memmax`   | memmax - the maximum amount of RAM (in GB) that terra can use when processing a raster dataset. Should be less than what is detected (see [`mem_info`](https://rspatial.github.io/terra/reference/mem.md), and higher values are ignored. Set it to a negative number or NA to ignore this value).                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `names`    | output layer names.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `NAflag`   | numeric. value to represent missing (`NA` or `NaN`) values. See note                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `scale`    | numeric. Cell values written to disk are divided by this value (default is 1). See [`scoff`](https://rspatial.github.io/terra/reference/scoff.md)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `offset`   | numeric. Value that is subtracted from the cell values written to disk (default is 0). See [`scoff`](https://rspatial.github.io/terra/reference/scoff.md)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `verbose`  | logical. If `TRUE` debugging information is printed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `steps`    | positive integers. In how many steps (chunks) do you want to process the data (for debugging)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `todisk`   | logical. If `TRUE` processing operates as if the dataset is very large and needs to be written to a temporary file (for debugging).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `metadata` | character, see `metags<-` to write metadata                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

## Note

GeoTiff files are, by default, written with LZW compression. If you do
not want compression, use `gdal="COMPRESS=NONE"`.

When writing integer values the lowest available value (given the
datatype) is used to represent `NA` for signed types, and the highest
value is used for unsigned values. This can be a problem with byte data
(between 0 and 255) as the value 255 is reserved for `NA`. To keep the
value 255, you need to set another value as `NAflag`, or do not set a
`NAflag` (with `NAflag=NA`)

## Examples

``` r
r <- rast(nrows=5, ncols=5, vals=1:25)

# create a temporary filename for the example
f <- file.path(tempdir(), "test.tif")

writeRaster(r, f, overwrite=TRUE)

writeRaster(r, f, overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"), datatype='INT1U')

## Or with a wopt argument:

writeRaster(r, f, overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='INT1U'))

## remove the file
unlink(f)
```
