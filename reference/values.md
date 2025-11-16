# Cell values and geometry attributes

Get the cell values of a SpatRaster or the attributes of a SpatVector.

By default all values returned are numeric. This is because a vector or
matrix can only store one data type, and a SpatRaster may consist of
multiple data types. However, if all layers have integer or logical
values, the returned values also have that datatype.

Note that with `values(x, dataframe=TRUE)` and
[`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md)`(x)`
the values returned match the type of each layer, and can be a mix of
numeric, logical, integer, and factor.

## Usage

``` r
# S4 method for class 'SpatRaster'
values(x, mat=TRUE, dataframe=FALSE, row=1, 
    nrows=nrow(x), col=1, ncols=ncol(x), na.rm=FALSE, ...)

# S4 method for class 'SpatVector'
values(x, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- mat:

  logical. If `TRUE`, values are returned as a matrix instead of as a
  vector, except when dataframe is `TRUE`

- dataframe:

  logical. If `TRUE`, values are returned as a `data.frame` instead of
  as a vector (also if matrix is `TRUE`)

- row:

  positive integer. Row number to start from, should be between 1 and
  nrow(x)

- nrows:

  positive integer. How many rows?

- col:

  positive integer. Column number to start from, should be between 1 and
  ncol(x)

- ncols:

  positive integer. How many columns? Default is the number of columns
  left after the start column

- na.rm:

  logical. Remove `NA`s?

- ...:

  additional arguments passed to
  [`data.frame`](https://rdrr.io/r/base/data.frame.html)

## Details

If `x` is a `SpatRaster`, and `mat=FALSE`, the values are returned as a
vector. In cell-order by layer. If `mat=TRUE`, a matrix is returned in
which the values of each layer are represented by a column (with
`ncell(x)` rows). The values per layer are in cell-order, that is, from
top-left, to top-right and then down by row. Use
[`as.matrix`](https://rspatial.github.io/terra/reference/coerce.md)`(x, wide=TRUE)`
for an alternative matrix representation where the number of rows and
columns matches that of `x`.

## Note

raster values that are `NA` (missing) are represented by `NaN`
(not-a-number) unless argument `dataframe` is `TRUE`.

## Value

matrix or data.frame

## See also

`values<-`,
[`focalValues`](https://rspatial.github.io/terra/reference/focalValues.md),
[`as.data.frame`](https://rspatial.github.io/terra/reference/as.data.frame.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : elev.tif 
#> name        : elevation 
#> min value   :       141 
#> max value   :       547 
x <- values(r)
x[3650:3655, ]
#> [1] 349 287 334 319 271 228
r[3650:3655]
#>   elevation
#> 1       349
#> 2       287
#> 3       334
#> 4       319
#> 5       271
#> 6       228


ff <- system.file("ex/lux.shp", package="terra")
v <- vect(ff)
y <- values(v)
head(y)
#>   ID_1       NAME_1 ID_2     NAME_2 AREA   POP
#> 1    1     Diekirch    1   Clervaux  312 18081
#> 2    1     Diekirch    2   Diekirch  218 32543
#> 3    1     Diekirch    3    Redange  259 18664
#> 4    1     Diekirch    4    Vianden   76  5163
#> 5    1     Diekirch    5      Wiltz  263 16735
#> 6    2 Grevenmacher    6 Echternach  188 18899
```
