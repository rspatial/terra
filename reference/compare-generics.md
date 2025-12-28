# Compare and logical methods

Standard comparison and logical operators for computations with
SpatRasters. Computations are local (applied on a cell by cell basis).
If multiple SpatRasters are used, these must have the same geometry
(extent and resolution). These operators have been implemented:

**Logical**: `!, &, |, isTRUE, isFALSE`

**Compare**:
` ==, !=, >, <, <=, >=, is.na, is.nan, is.finite, is.infinite`

See [`not.na`](https://rspatial.github.io/terra/reference/not.na.md) for
the inverse of `is.na`, and
[`noNA`](https://rspatial.github.io/terra/reference/summarize-generics.md)
to detect cells with missing value across layers.

The `compare` and `logic` methods implement these operators in a method
that can return `NA` istead of `FALSE` and allows for setting an output
filename.

The terra package does not distinguish between `NA` (not available) and
`NaN` (not a number). In most cases this state is represented by `NaN`.

If you use a SpatRaster with a vector of multiple numbers, each element
in the vector is considered a layer (with a constant value). If you use
a SpatRaster with a matrix, the number of columns of the matrix must
match the number of layers of the SpatRaster. The rows are used to match
the cells. That is, if there are two rows, these match cells 1 and 2,
and they are recycled to 3 and 4, etc.

The following method has been implemented for **(SpatExtent,
SpatExtent)**: `==`

## Usage

``` r
# S4 method for class 'SpatRaster'
compare(x, y, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatRaster'
logic(x, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...)
```

## Arguments

- x:

  SpatRaster

- y:

  SpatRaster or numeric

- oper:

  character. Operator name. For `compare` this can be one of
  `"==", "!=", ">", "<", ">=", "<="` and for `logic` it can be one of
  `"!", "is.na", "not.na", "allNA", "anyNA", "noneNA", "is.infinite", "is.finite", "iSTRUE", "isFALSE"`

- falseNA:

  logical. Should the result be `TRUE, NA` instead of `TRUE, FALSE`?

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`all.equal`](https://rspatial.github.io/terra/reference/all.equal.md),
[Arith-methods](https://rspatial.github.io/terra/reference/arith-generic.md).
See [`ifel`](https://rspatial.github.io/terra/reference/ifelse.md) to
conveniently combine operations and
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md)
or [`app`](https://rspatial.github.io/terra/reference/app.md) to apply
any R function to a SpatRaster.

## Value

SpatRaster or SpatExtent

## Examples

``` r
r1 <- rast(ncols=10, nrows=10)
values(r1) <- runif(ncell(r1))
r1[10:20] <- NA
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)

x <- is.na(r1)
!x
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   : FALSE 
#> max value   :  TRUE 
r1 == r2 
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   : FALSE 
#> max value   : FALSE 
compare(r1, r2, "==")
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   : FALSE 
#> max value   : FALSE 
compare(r1, r2, "==", TRUE)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        : lyr.1 
#> min value   :  TRUE 
#> max value   :  TRUE 
```
