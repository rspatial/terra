# Data type of a SpatRaster or SpatVector

Get the data types of the fields (attributes, variables) of a SpatVector
or of the file(s) associated with a SpatRaster. A (layer of a)
SpatRaster has no datatype if it has no values, or if the values are in
memory.

## Usage

``` r
# S4 method for class 'SpatRaster'
datatype(x, bylyr=TRUE)

# S4 method for class 'SpatVector'
datatype(x)
```

## Arguments

- x:

  SpatRaster or SpatVector

- bylyr:

  logical. If `TRUE` a value is returned for each layer. Otherwise, a
  value is returned for each data source (such as a file)

## Details

Setting the data type is useful if you want to write values to disk with
[`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md).
In other cases you can use functions such as `round` and `floor`, or
[`as.bool`](https://rspatial.github.io/terra/reference/is.bool.md)

raster datatypes are described by 5 characters. The first three indicate
whether the values are integer or decimal values. The fourth character
indicates the number of bytes used to save the values on disk, and the
last character indicates whether the numbers are signed (that is, can be
negative and positive values) or not (only zero and positive values
allowed)

The following raster datatypes are available:

|                         |                            |                            |
|-------------------------|----------------------------|----------------------------|
| **Datatype definition** | **minimum possible value** | **maximum possible value** |
| `INT1U`                 | 0                          | 255                        |
| `INT2U`                 | 0                          | 65,535                     |
| `INT4U`                 | 0                          | 4,294,967,296              |
| `INT8U`                 | 0                          | 18,446,744,073,709,551,616 |
| `INT1S`                 | -128                       | 128                        |
| `INT2S`                 | -32,767                    | 32,767                     |
| `INT4S`                 | -2,147,483,647             | 2,147,483,647              |
| `INT8S`                 | -9,223,372,036,854,775,808 | 9,223,372,036,854,775,808  |
| `FLT4S`                 | -3.4e+38                   | 3.4e+38                    |
| `FLT8S`                 | -1.7e+308                  | 1.7e+308                   |

For all integer and byte types the lowest (signed) or highest (unsigned)
value is used to store `NA`. For float types `NaN` is used (following
the IEEE 754 standard).

Note that very large integer numbers may be imprecise as they are
internally represented as decimal numbers.

Also note that `NaN` may not be equally supported by all
implementations. For example OGR SQL and SQLite queries generally
convert `NaN` values to `NULL`.

`INT4U` and `INT8U` are available but they are best avoided as R does
not support 32-bit or 64-bit unsigned integers. `INT8U` is a special
case where the NA store value is `18446744073709549568`
(`UINT64_MAX - 1101`) because of precision in decimal representation.

`INT8U` and `INT8S` require GDAL version 3.5 or higher. `INT1S` requires
GDAL version 3.7 or higher.

## Value

character

## See also

[`Raster data types`](https://rspatial.github.io/terra/reference/is.bool.md)
to check / set the type of SpatRaster values.

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
datatype(v)
#> [1] "double" "string" "double" "string" "double" "double"

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
datatype(r)
#> [1] "INT2S"

# no data type
datatype(rast()) 
#> [1] ""
```
