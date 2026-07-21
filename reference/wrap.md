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
#> Error in paste0("(", x, ")"): non-string argument to .Internal(paste0)
p
#> Error: object 'p' not found
vv <- vect(p)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'vect': object 'p' not found
vv
#> Error: object 'vv' not found
```
