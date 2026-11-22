# Subdivide a raster or polygons

Divide a `SpatRaster` by splitting it `n` times. By default the
resulting parts have approximately the same sum of cell values
(`split="values"`). Use `split="notNA"` for the same number of non-`NA`
cells, or `split="cells"` to ignore values and `NA`s and split the grid
by cell count.

Divide a `SpatVector` of polygons into `n` compact and approximately
equal-area parts. The results are not deterministic, so you should use
[`set.seed`](https://rdrr.io/r/base/Random.html) to be able to reproduce
your results. If you get a warning about non-convergence, you can
increase the number of iterations used with additional argument
`iter.max`

## Usage

``` r
# S4 method for class 'SpatRaster'
divide(x, n=2, start="ns", as.raster=FALSE, na.rm=TRUE, split="values")

# S4 method for class 'SpatVector'
divide(x, n=5, w=NULL, alpha=1, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector of polygons

- n:

  numeric.

  For a SpatVector, a single positive integer: the number of parts.

  For a SpatRaster, a single positive integer is the number of split
  rounds. Each round cuts every current piece in two, alternating
  direction (see `start`). So `n=1` gives 2 parts, `n=2` gives 4, `n=3`
  gives 8.

  For a SpatRaster, `n` can instead be a vector that lists the cuts, in
  order. Each value is applied to every piece from the previous step.
  The sign is the direction: positive = North-South (a horizontal cut),
  negative = East-West (a vertical cut). The absolute value is the
  number of cuts in that step, so the piece is divided into `abs(n)+1`
  parts (`1` = halves, `2` = thirds, `3` = quarters, ...). For example
  `c(1, -1)` is a North-South cut into 2, then an East-West cut of each
  half (4 parts); `c(2, -1, 1)` is a North-South cut into 3, then
  East-West into 2, then North-South into 2 (12 parts).

  For a SpatVector, `n` can be a list with at least one of these
  elements: `horizontal` and `vertical` that specify the proportions of
  the area that splits should cover. This can either be a single
  fraction such as 1/3, or a sequence of fractions in ascending order
  such as `c(1/4, 1/2, 1)`

- start:

  character. Initial split direction when `n` is a single number. `"ns"`
  for North-South (horizontal cuts) or `"ew"` for East-West (vertical
  cuts). Ignored if `n` is a vector

- as.raster:

  logical. If `FALSE` a SpatVector is returned. If `TRUE`, a SpatRaster
  is returned. If `NA` a list with a SpatRaster and a SpatVector is
  returned

- na.rm:

  logical. If `TRUE` cells in `x` that are `NA` are not included in the
  output

- split:

  character. For SpatRaster only. `"values"` (the default) balances the
  sum of cell values in each part. `"notNA"` balances the number of
  non-`NA` cells. `"cells"` ignores values and `NA`s and splits by the
  number of cells

- w:

  SpatRaster with, for example, environmental data

- alpha:

  numeric. One or two numbers that act as weights for the x and y
  coordinates

- ...:

  additional arguments such as `iter.max` passed on to
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## See also

[`thresh`](https://rspatial.github.io/terra/reference/thresh.md);
[`makeTiles`](https://rspatial.github.io/terra/reference/makeTiles.md);
[`split`](https://rspatial.github.io/terra/reference/split.md);
[`segregate`](https://rspatial.github.io/terra/reference/segregate.md)

## Value

SpatVector or SpatRaster, or a list with both

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- divide(r, 3)
# plot(r); lines(x)


f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
d <- divide(v, 3)
dv <- divide(v, list(h=.5))
```
