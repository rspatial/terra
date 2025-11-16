# wrap and unwrap

Use `wrap` to pack a SpatVector or SpatRaster\* to create a Packed\*
object. Packed objects can be passed over a connection that serializes
(e.g. to nodes on a computer cluster). At the receiving end they need to
be unpacked with `unwrap`.

## Usage

``` r
# S4 method for class 'SpatRaster'
wrap(x, proxy=FALSE)

# S4 method for class 'SpatRasterDataset'
wrap(x, proxy=FALSE)

# S4 method for class 'SpatRasterCollection'
wrap(x, proxy=FALSE)

# S4 method for class 'SpatVector'
wrap(x)

# S4 method for class 'ANY'
unwrap(x)
```

## Arguments

- x:

  SpatVector, SpatRaster, SpatRasterDataset or SpatRasterCollection

- proxy:

  logical. If `FALSE` raster cell values are forced to memory if
  possible. If `TRUE`, a reference to source filenames is stored for
  data sources that are not in memory

## Value

`wrap`: Packed\* object

`unwrap`: SpatVector, SpatRaster, SpatRasterCollection,
SpatRasterDataset

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
p <- wrap(v)
p
#> [1] "This is a PackedSpatVector object. Use 'terra::unwrap()' to unpack it"
vv <- vect(p)
vv
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
