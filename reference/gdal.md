# GDAL version, supported file formats, and cache size

Set the `GDAL` warning level or get a `data.frame` with the available
GDAL drivers (file formats), or, if `warn=NA` and `drivers=FALSE`, you
get the version numbers of one or all of the GDAL, PROJ and GEOS
libraries.

`GDAL` is the software library that terra builds on to read and write
spatial data and for some raster data processing. `PROJ` is used for
transformation of coordinates ("projection") and `GEOS` is used for
geometric operations with vector data.

The current `GDAL` configuration options and obtained with
`getGDALconfig` and changed with `setGDALconfig`.

## Usage

``` r
gdal(warn=NA, drivers=FALSE, ...)
gdalCache(size=NA)
setGDALconfig(option, value="")
getGDALconfig(option)
clearVSIcache()
libVersion(lib="all", parse=FALSE)
unloadGDALdrivers(x)
```

## Arguments

- warn:

  If `NA` and `drivers=FALSE`, the version of the library specified by
  `lib` is returned. Otherwise, the value should be an integer between 1
  and 4 representing the level of GDAL warnings and errors that are
  passed to R. 1 = warnings and errors; 2 = errors only (recoverable
  errors as a warning); 3 = irrecoverable errors only; 4 = ignore all
  errors and warnings. The default setting is 2

- drivers:

  logical. If `TRUE` a data.frame with the raster and vector data
  formats that are available.

- ...:

  additional arguments (for backwards compatibility only)

- size:

  numeric. The new cache size in MB

- option:

  character. GDAL configuration option name, or a "name=value" string
  (in which case the value argument is ignored

- value:

  character. value for GDAL configuration option. Use "" to reset it to
  its default value

- lib:

  character. "gdal", "proj", or "geos", or any other value to get the
  versions numbers of all three

- parse:

  logical. Should the version be parsed into three numerical values
  (major, minor and sub versions)?

- x:

  character. Drivers names such as "GTiff" to be unloaded. Or "" to
  reload all drivers

## See also

[`describe`](https://rspatial.github.io/terra/reference/describe.md) for
file-level metadata "GDALinfo"

## Value

character

## Examples

``` r
gdal()
#> [1] "3.8.4"
gdal(2)
head(gdal(drivers=TRUE))
#>      name raster vector        can  vsi                     long.name
#> 1 AAIGrid   TRUE  FALSE read/write TRUE           Arc/Info ASCII Grid
#> 2    ACE2   TRUE  FALSE       read TRUE                          ACE2
#> 3    ADRG   TRUE  FALSE read/write TRUE ARC Digitized Raster Graphics
#> 4     AIG   TRUE  FALSE       read TRUE          Arc/Info Binary Grid
#> 5     ARG   TRUE  FALSE read/write TRUE     Azavea Raster Grid format
#> 6  AVCBin  FALSE   TRUE       read TRUE      Arc/Info Binary Coverage
libVersion("all", TRUE)
#>      major minor sub
#> gdal     3     8   4
#> proj     9     4   0
#> geos     3    12   1
```
