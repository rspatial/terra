# Update a raster file

Update metadata or cell values of a file that is the data source of a
SpatRaster. This modifies the file on disk directly without reading the
entire raster into memory.

BE CAREFUL with this function as you are overwriting data in an existing
file.

## Usage

``` r
# S4 method for class 'SpatRaster'
update(object, crs=FALSE, extent=FALSE, names=FALSE,
    cells=NULL, values=NULL, layer=0)
```

## Arguments

- object:

  SpatRaster with a file source

- crs:

  logical. Should the coordinate reference system be updated?

- extent:

  logical. Should the extent be updated?

- names:

  logical. Should the names be updated?

- cells:

  numeric. Cell numbers to update (1-indexed). Must be used together
  with `values`

- values:

  numeric. New values for the specified cells. Can be a single value
  (recycled), one value per cell (applied to all layers), or
  `length(cells) * length(layer)` values (one per cell per layer, in
  layer-major order)

- layer:

  numeric or character. Layer(s) to update. Use `0` (default) to update
  all layers. Positive integers or layer names select specific layers

## Value

SpatRaster (invisibly)

## Details

The `cells` and `values` arguments allow updating specific cell values
on disk without loading the raster into memory. This is useful for
modifying very large rasters that cannot fit in memory. See
[`set.values`](https://rspatial.github.io/terra/reference/inplace.md)
for modifying values of in-memory rasters.

## Examples

``` r
## update metadata
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
#> source      : file21af336e791e.tif
#> names       :   A,   B,   C
#> min values  :   0,   0,   0
#> max values  : 255, 255, 255

## update cell values on disk
r <- rast(nrows=10, ncols=10, vals=1:100)
f <- paste0(tempfile(), ".tif")
x <- writeRaster(r, f)
x[c(1, 2, 50, 51, 100)]
#>   lyr.1
#> 1     1
#> 2     2
#> 3    50
#> 4    51
#> 5   100
update(x, cells=c(1, 50, 100), values=c(999, 888, 777))
rast(f)[c(1, 2, 50, 51, 100)]
#>   lyr.1
#> 1   999
#> 2     2
#> 3   888
#> 4    51
#> 5   777

## update specific layer only
r <- rast(nrows=5, ncols=5, nlyr=3, vals=1:75)
f <- paste0(tempfile(), ".tif")
x <- writeRaster(r, f, datatype="FLT4S")
x[c(1, 2, 25, 26)]
#>   lyr.1 lyr.2 lyr.3
#> 1     1    26    51
#> 2     2    27    52
#> 3    25    50    75
#> 4    NA    NA    NA
update(x, cells=c(1, 25), values=c(0, 0), layer=2)
rast(f)[c(1, 2, 25, 26)]
#>   lyr.1 lyr.2 lyr.3
#> 1     1     0    51
#> 2     2    27    52
#> 3    25     0    75
#> 4    NA    NA    NA
```
