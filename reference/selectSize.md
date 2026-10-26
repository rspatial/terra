# Select features by size

Select the largest or smallest polygons or lines in a `SpatVector`,
optionally providing a maximum number of selected features, filtering by
a size threshold, dissolving adjacent features, and/or splitting
multipart features into single parts.

## Usage

``` r
selectSize(x, n=Inf, largest=TRUE, threshold=NULL,
           dissolve=FALSE, singleparts=TRUE)
```

## Arguments

- x:

  `SpatVector` of polygons or lines

- n:

  optional integer value (default Inf) indicating the maximum number of
  largest or smallest features to select (if they pass the optional
  threshold)

- largest:

  if `TRUE` (default) select the largest feature(s), otherwise the
  smallest

- threshold:

  optional size threshold (area in square meters for polygons, length in
  meters for lines) above which (if `largest=TRUE` or below which (if
  `largest=FALSE` features should be selected

- dissolve:

  if `TRUE`, aggregate features before selection. Default `FALSE`

- singleparts:

  if `TRUE` (default), split multipart features into individual parts
  before selection

## Details

Polygon size is calculated using
[`expanse`](https://rspatial.github.io/terra/reference/expanse.md); line
size is calculated using
[`perim`](https://rspatial.github.io/terra/reference/perim.md). The
threshold is applied before selecting the requested number of features.
If fewer than `n` features satisfy the threshold, only qualifying
features are returned.

If `dissolve=TRUE`, features are aggregated before size calculation. If
`singleparts=TRUE`, multipart features are disaggregated before size
calculation.

## Value

`SpatVector` containing the selected features.

## Author

A. Márcia Barbosa

## See also

[`expanse`](https://rspatial.github.io/terra/reference/expanse.md),
[`perim`](https://rspatial.github.io/terra/reference/perim.md),
[`selectHighest`](https://rspatial.github.io/terra/reference/selectHigh.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

# Select the five largest polygons
large5 <- selectSize(v, n=5)

# Select the three smallest polygons
small3 <- selectSize(v, n=3, largest=FALSE)

# Select up to 4 largest polygons satisfying a given minimum area
large <- selectSize(v, n=4, threshold=260e+6)

if (FALSE) { # \dontrun{
plot(v)
polys(large5, col="red")
polys(small3, col="blue")
polys(large, col="gold", density=50)
} # }
```
