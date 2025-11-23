# Options

Get or set general options.

## Usage

``` r
terraOptions(..., print=TRUE)
```

## Arguments

- ...:

  option names and values (see Details). Or missing, to get or show the
  current options

- print:

  logical. If `TRUE` the option names and values are printed

## Details

The following options are available.

**memfrac** - value between 0 and 0.9 (larger values give a warning).
The fraction of RAM that may be used by the program.

**memmin** - if memory required is below this threshold (in GB), the
memory is assumed to be available. Otherwise, terra checks if it is
available.

**memmax** - the maximum amount of RAM (in GB) that terra is allowed to
use when processing a raster dataset. Should be less than what is
detected (see
[`mem_info`](https://rspatial.github.io/terra/reference/mem.md)), and
higher values are ignored. Set it to a negative number or NA to not set
this option. `terraOptions` only shows the value of `memmax` if it is
set.

**tempdir** - directory where temporary files are written. The default
what is returned by [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

**datatype** - default data type. See
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md).

**todisk** - logical. If `TRUE` write all raster data to disk (temp file
if no file name is specified). For debugging.

**progress** - non-negative integer. A progress bar is shown if the
number of chunks in which the data is processed is larger than this
number. No progress bar is shown if the value is zero.

**verbose** - logical. If `TRUE` debugging info is printed for some
functions.

**tolerance** - numeric. Difference in raster extent (expressed as the
fraction of the raster resolution) that can be ignored when comparing
alignment of rasters.

## Note

It is possible to set your own default options in "etc/.Rprofile.site"
of your R installation like this

`options(terra_default=list(tempdir="d:/temp", memfrac=.4))`

But that may not be a good practice. It is clearer to set your favorite
options at the beginning of each script.

## Value

list. Invisibly if `print=TRUE`

## Examples

``` r
terraOptions()
#> memfrac   : 0.5
#> tolerance : 0.1
#> verbose   : FALSE
#> todisk    : FALSE
#> tempdir   : /tmp/Rtmps0DlvB
#> datatype  : FLT4S
#> memmin    : 1
#> progress  : 3
terraOptions(memfrac=0.5, tempdir = "c:/temp")
#> Warning: [options] you cannot set the tempdir to a path that does not exist
terraOptions(progress=10)
terraOptions()
#> memfrac   : 0.5
#> tolerance : 0.1
#> verbose   : FALSE
#> todisk    : FALSE
#> tempdir   : /tmp/Rtmps0DlvB
#> datatype  : FLT4S
#> memmin    : 1
#> progress  : 10
```
