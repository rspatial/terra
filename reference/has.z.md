# Does a SpatVector have Z coordinates?

Returns `TRUE` if any geometry part (or polygon hole) stores Z values.

Z storage is optional: empty Z means a conventional 2D geometry. When Z
is present it is included in
[`geom`](https://rspatial.github.io/terra/reference/geometry.md) and
[`crds`](https://rspatial.github.io/terra/reference/crds.md), preserved
by
[`writeVector`](https://rspatial.github.io/terra/reference/writeVector.md)
/ [`vect`](https://rspatial.github.io/terra/reference/vect.md) file I/O
(25D geometries), and by coercion to/from sf (see also
<https://github.com/rspatial/terra/issues/824>).

## Usage

``` r
# S4 method for class 'SpatVector'
has.z(x)
```

## Arguments

- x:

  SpatVector

## Value

logical

## See also

[`geom`](https://rspatial.github.io/terra/reference/geometry.md)`, `[`crds`](https://rspatial.github.io/terra/reference/crds.md)`, `[`vect`](https://rspatial.github.io/terra/reference/vect.md)

## Examples

``` r
v <- vect(cbind(1:3, 1:3))
has.z(v)
#> [1] FALSE
```
