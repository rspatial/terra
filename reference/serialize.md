# saveRDS and serialize for SpatVector and SpatRaster\*

serialize and saveRDS for SpatVector, SpatRaster, SpatRasterDataset and
SpatRasterCollection. Note that these objects will first be "packed"
with [`wrap`](https://rspatial.github.io/terra/reference/wrap.md), and
after unserialize/readRDS they need to be unpacked with `rast` or
`vect`.

Extensive use of these functions is not recommended. Especially for
SpatRaster it is generally much more efficient to use
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)
and write, e.g., a GTiff file.

## Usage

``` r
# S4 method for class 'SpatRaster'
saveRDS(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

# S4 method for class 'SpatRasterDataset'
saveRDS(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

# S4 method for class 'SpatRasterCollection'
saveRDS(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

# S4 method for class 'SpatVector'
saveRDS(object, file="", ascii = FALSE, version = NULL, compress=TRUE, refhook = NULL)

# S4 method for class 'SpatRaster'
serialize(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)

# S4 method for class 'SpatVector'
serialize(object, connection, ascii = FALSE, xdr = TRUE, version = NULL, refhook = NULL)
```

## Arguments

- object:

  SpatVector, SpatRaster, SpatRasterDataset or SpatRasterCollection

- file:

  file name to save object to

- connection:

  see `serialize`

- ascii:

  see `serialize` or `saveRDS`

- version:

  see `serialize` or `saveRDS`

- compress:

  see `serialize` or `saveRDS`

- refhook:

  see `serialize` or `saveRDS`

- xdr:

  see `serialize` or `saveRDS`

## Value

Packed\* object

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
p <- serialize(v, NULL)
head(p)
#> [1] 58 0a 00 00 00 03
x <- unserialize(p)
x
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 12, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.44781, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04
```
