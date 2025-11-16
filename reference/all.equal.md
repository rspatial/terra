# Compare two SpatRaster, SpatVector, or SpatExtent objects for equality

Compare two objects for (near) equality

In the case of SpatRasters, first the attributes of the objects are
compared. If these are the same, a (perhaps small) sample of the raster
cells is compared as well.

The sample size used can be increased with the `maxcell` argument. You
can set it to `Inf`, but for large rasters your computer may not have
sufficient memory. See the examples for a safe way to compare all
values.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
all.equal(target, current, maxcell=100000, ...)

# S4 method for class 'SpatVector,SpatVector'
all.equal(target, current, ...)

# S4 method for class 'SpatExtent,SpatExtent'
all.equal(target, current, ...)
```

## Arguments

- target:

  SpatRaster, SpatVector, or SpatExtent

- current:

  object of the same class as `target`

- maxcell:

  positive integer. The size of the regular sample used to compare cell
  values

- ...:

  additional arguments passed to
  [`all.equal.numeric`](https://rdrr.io/r/base/all.equal.html) to
  compare cell values for SpatRaster and geometry and attribute values
  for SpatVectors

## See also

[`identical`](https://rspatial.github.io/terra/reference/identical.md),
[`compareGeom`](https://rspatial.github.io/terra/reference/compareGeom.md)

## Value

Either `TRUE` or a character vector describing the differences between
target and current.

## Examples

``` r
x <- sqrt(1:100)
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = mat)

all.equal(r1, r2)
#> [1] "Mean relative difference: 0.3858482"
all.equal(r1, r1*1)
#> [1] TRUE
all.equal(rast(r1), rast(r2))
#> [1] TRUE

# compare geometries 
compareGeom(r1, r2)
#> [1] TRUE

# Compare all cell values for near equality
# as floating point number imprecision can be a problem
m <- minmax(r1 - r2)
all(abs(m) < 1e-7)
#> [1] FALSE

# comparison of cell values to create new SpatRaster
e <- r1 == r2
```
