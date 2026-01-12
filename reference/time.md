# time of SpatRaster layers

Get or set the time of the layers of a SpatRaster. Time can be stored as
[`POSIXlt`](https://rdrr.io/r/base/DateTimeClasses.html) (date and time,
with a resolution of seconds, and a time zone),
[`Date`](https://rdrr.io/r/base/Dates.html), "months", "years", or
"yearmonths".

`timeInfo` and `has.time` are helper functions to understand what a time
data a SpatRaster has.

## Usage

``` r
# S4 method for class 'SpatRaster'
has.time(x)

# S4 method for class 'SpatRaster'
time(x, format="")

# S4 method for class 'SpatRaster'
time(x, tstep = "") <- value

# S4 method for class 'SpatRaster'
timeInfo(x)
```

## See also

[`depth`](https://rspatial.github.io/terra/reference/depth.md)

## Arguments

- x:

  SpatRaster or SpatRasterDataset

- format:

  One of "", "seconds" (POSIXlt), "days" (Date), "yearmonths" (decimal
  years), "years", "months". If "", the returned format is (based on)
  the format that was used to set the time

- value:

  `Date`, `POSIXt`, `yearmon` (defined in package zoo), or numeric

- tstep:

  One of "years", "months", "yearmonths". Used when `value` is numeric.
  Ignored when `value` is of type `Date`, `POSIXt`, or `yearmon`

## Value

`time`: POSIXlt, Date, or numeric `timeInfo`: data.frame with time step
and time zone information (if available) `has.time`: logical

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   

# Date"
d <- as.Date("2001-05-04") + 0:2
time(s) <- d
time(s)
#> [1] "2001-05-04" "2001-05-05" "2001-05-06"

# POSIX (date/time with a resolution of seconds)
time(s) <- as.POSIXlt(d)
time(s)
#> [1] "2001-05-04 UTC" "2001-05-05 UTC" "2001-05-06 UTC"

# with time zone
time(s) <- as.POSIXlt(Sys.time(), "America/New_York") + 0:2
time(s)
#> [1] "2026-01-12 01:30:36 EST" "2026-01-12 01:30:37 EST"
#> [3] "2026-01-12 01:30:38 EST"
timeInfo(s)
#>   time    step             zone
#> 1 TRUE seconds America/New_York

# years
time(s, tstep="years") <- 2000 + 0:2
s
#> class       : SpatRaster 
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> colors RGB  : 1, 2, 3 
#> names       : red, green, blue 
#> min values  :   0,     0,    0 
#> max values  : 255,   255,  255 
#> time (years): 2000 to 2002 (3 steps) 

time(s, tstep="months") <- 1:3
s 
#> class       : SpatRaster 
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> colors RGB  : 1, 2, 3 
#> names       : red, green, blue 
#> min values  :   0,     0,    0 
#> max values  : 255,   255,  255 
#> time (mnts) : Jan to Mar (3 steps) 
```
