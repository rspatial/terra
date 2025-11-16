# Arithmetic

Standard arithmetic operators for computations with SpatRasters.
Computations are local (applied on a cell by cell basis). If multiple
SpatRasters are used, these must have the same geometry (extent and
resolution). These operators have been implemented:

` +, -, *, /, ^, %%, %/% `

You can also use a SpatRaster and a vector or a matrix. If you use a
SpatRaster with a vector of multiple numbers, each element in the vector
is considered a layer (with a constant value). If you use a SpatRaster
with a matrix, the number of columns of the matrix must match the number
of layers of the SpatRaster. The rows are used to match the cells. That
is, if there are two rows, these match cells 1 and 2, and they are
recycled to 3 and 4, etc.

The following methods have been implemented for
`(SpatExtent, SpatExtent)`: `+, -`, and the following for
`(SpatExtent, numeric)`: `+, -, *, /, %%`

## See also

[`ifel`](https://rspatial.github.io/terra/reference/ifelse.md) to
conveniently combine operations and
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md)
or [`app`](https://rspatial.github.io/terra/reference/app.md) to use
mathematical functions not implemented by the package.

## Value

SpatRaster or SpatExtent

## Examples

``` r
r1 <- rast(ncols=10, nrows=10)
v <- runif(ncell(r1))
v[10:20] <- NA
values(r1) <- v
r2 <- rast(r1)
values(r2) <- 1:ncell(r2) / ncell(r2)
r3 <- r1 + r2
r2 <- r1 / 10
r3 <- r1 * (r2 - 1 / r2)

b <- c(r1, r2, r3)
b2 <- b * 10

### SpatExtent methods
x <- ext(0.1, 2.2, 0, 3)
y <- ext(-2, 1, -2,2)
# union
x + y
#> SpatExtent : -2, 2.2, -2, 3 (xmin, xmax, ymin, ymax)
# intersection
x * y
#> SpatExtent : 0.1, 1, 0, 2 (xmin, xmax, ymin, ymax)

e <- x %% 2
e
#> SpatExtent : 0, 4, 0, 4 (xmin, xmax, ymin, ymax)
e * 2
#> SpatExtent : -2, 6, -2, 6 (xmin, xmax, ymin, ymax)
e / 2
#> SpatExtent : 1, 3, 1, 3 (xmin, xmax, ymin, ymax)
e + 1
#> SpatExtent : -1, 5, -1, 5 (xmin, xmax, ymin, ymax)
e - 1
#> SpatExtent : 1, 3, 1, 3 (xmin, xmax, ymin, ymax)
```
