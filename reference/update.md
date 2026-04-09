# Change values in a file

Change the contents of a file that is the data source of a SpatRaster.
BE CAREFUL as you are overwriting values in an existing file.

## Usage

``` r
# S4 method for class 'SpatRaster'
update(object, crs=FALSE, extent=FALSE, names=FALSE)
```

## Arguments

- object:

  SpatRaster

- crs:

  logical. Should the coordinate reference system be updated?

- extent:

  logical. Should the extent be updated?

- names:

  logical. Should the names be updated?

## Value

SpatRaster (invisibly)

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s*1, fname)

ext(x) <- ext(x) + 1
crs(x) <- "+proj=utm +zone=1"
names(x) <- LETTERS[1:3]

update(x, crs=TRUE, extent=TRUE, names=TRUE)

rast(fname)
#> class       : SpatRaster
#> size        : 77, 101, 3  (nrow, ncol, nlyr)
#> resolution  : 1.019802, 1.025974  (x, y)
#> extent      : -1, 102, -1, 78  (xmin, xmax, ymin, ymax)
#> coord. ref. : WGS 84 / UTM zone 1N (EPSG:32601)
#> source      : file21625aecc03b.tif
#> names       :   A,   B,   C
#> min values  :   0,   0,   0
#> max values  : 255, 255, 255
```
