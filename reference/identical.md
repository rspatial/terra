# Compare two SpatRaster, SpatVector or SpatExtent objects for equality

When, comparing two SpatRasters for equality, first the attributes of
the objects are compared. If these are the same, a the raster cells are
compared as well. This can be time consuming, and you may prefer to use
a sample instead with
[`all.equal`](https://rspatial.github.io/terra/reference/all.equal.md)

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
identical(x, y)

# S4 method for class 'SpatVector,SpatVector'
identical(x, y)

# S4 method for class 'SpatExtent,SpatExtent'
identical(x, y)
```

## Arguments

- x:

  SpatRaster, SpatVector, or SpatExtent

- y:

  object of the same class as `x`

## See also

[`all.equal`](https://rspatial.github.io/terra/reference/all.equal.md),
[`compareGeom`](https://rspatial.github.io/terra/reference/compareGeom.md)

## Value

single logical value

## Examples

``` r
x <- sqrt(1:100)
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = t(mat))

identical(r1, r2)
#> [1] TRUE
identical(r1, r1*1)
#> [1] TRUE
identical(rast(r1), rast(r2))
#> [1] TRUE
```
