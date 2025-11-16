# Write SpatVector data to a file

Write a SpatVector to a file. You can choose one of many file formats.

## Usage

``` r
# S4 method for class 'SpatVector,character'
writeVector(x, filename, filetype=NULL, layer=NULL, insert=FALSE,
    overwrite=FALSE, options="ENCODING=UTF-8")
```

## Arguments

- x:

  SpatVector

- filename:

  character. Output filename

- filetype:

  character. A file format associated with a GDAL "driver" such as "ESRI
  Shapefile". See `gdal(drivers=TRUE)` or the [GDAL
  docs](https://gdal.org/en/latest/drivers/vector/index.html). If `NULL`
  it is attempted to guess the filetype from the filename extension

- layer:

  character. Output layer name. If `NULL` the filename is used

- insert:

  logical. If `TRUE`, a new layer is inserted into the file, or an
  existing layer overwritten (if `overwrite=TRUE`), if the format
  supports it (e.g. GPKG allows that). See
  [`vector_layers`](https://rspatial.github.io/terra/reference/vector_layers.md)
  to remove a layer

- overwrite:

  logical. If `TRUE` and `insert=FALSE`, `filename` is overwritten if
  the file format and layer structure permits it. If `TRUE` and
  `insert=TRUE`, only the target layer is overwritten if the format
  supports it (e.g. GPKG).

- options:

  character. Format specific GDAL options such as "ENCODING=UTF-8". Use
  NULL or "" to not use any options

## Examples

``` r
v <- vect(cbind(1:5,1:5))
crs(v) <- "+proj=longlat +datum=WGS84"
v$id <- 1:length(v)
v$name <- letters[1:length(v)]
tmpf1 <- paste0(tempfile(), ".gpkg")
writeVector(v, tmpf1, overwrite=TRUE)
x <- vect(tmpf1)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
tmpf2 <- paste0(tempfile(), ".gpkg")
writeVector(v, tmpf2, overwrite=TRUE)
y <- vect(tmpf2)
```
