# Raster value types

The values in a SpatRaster layer are by default numeric, but they can
also be set to be logical (Boolean), integer, or categorical (factor).

For a `SpatRaster`, `as.logical` and `isTRUE` is equivalent to
`as.bool`. `isFALSE` is equivalent to `!as.bool`, and `as.integer` is
the same as `as.int`.

`as.bool` and `as.int` force the values into the correct range (e.g.
whole integers) but in-memory cell values are still stored as numeric.
They will behave like the assigned types, though, and will be written to
files with that data type (if the file type supports it).

See [`levels`](https://rspatial.github.io/terra/reference/factors.md)
and [`cats`](https://rspatial.github.io/terra/reference/factors.md) to
create categorical layers by setting labels.

## Usage

``` r
# S4 method for class 'SpatRaster'
is.num(x)

# S4 method for class 'SpatRaster'
is.bool(x)

# S4 method for class 'SpatRaster'
as.bool(x, filename, ...)

# S4 method for class 'SpatRaster'
is.int(x)

# S4 method for class 'SpatRaster'
as.int(x, filename, ...)

# S4 method for class 'SpatRaster'
is.factor(x)

# S4 method for class 'SpatRaster'
as.factor(x)
```

## Arguments

- x:

  SpatRaster

- filename:

  character. Output filename

- ...:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`levels`](https://rspatial.github.io/terra/reference/factors.md) and
[`cats`](https://rspatial.github.io/terra/reference/factors.md) to
create categorical layers (and set labels).

## Value

The `as.*` methods return a new `SpatRaster`, whereas the `is.*` methods
return a `logical` value for each layer in `x`.

## Examples

``` r
r <- rast(nrows=10, ncols=10, vals=1:100)
is.bool(r)
#> [1] FALSE
z <- as.bool(r)
is.bool(z)
#> [1] TRUE

x <- r > 25
is.bool(x)
#> [1] TRUE

rr <- r/2
is.int(rr)
#> [1] FALSE
is.int(round(rr))
#> [1] TRUE
```
