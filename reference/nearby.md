# nearby geometries

Identify geometries that are near to each other. Either get the index of
all geometries within a certain distance, or the k nearest neighbors, or
(with `nearest`) get the nearest points between two geometries.

## Usage

``` r
# S4 method for class 'SpatVector'
nearby(x, y=NULL, distance=0, k=1, centroids=TRUE, symmetrical=TRUE, method="geo")

# S4 method for class 'SpatVector'
nearest(x, y, pairs=FALSE, centroids=TRUE, lines=FALSE, method="geo")
```

## Arguments

- x:

  SpatVector

- y:

  SpatVector or NULL

- distance:

  numeric. maximum distance

- k:

  positive integer. number of neighbors. Ignored if `distance > 0`

- centroids:

  logical. Should the centroids of polygons be used?

- symmetrical:

  logical. If `TRUE`, a near pair is only included once. That is, if
  geometry 1 is near to geometry 3, the implied nearness between 3 and 1
  is not reported. Ignored if `k` neighbors are returned

- method:

  character. One of "geo", "haversine", "cosine". With "geo" the most
  precise but slower method of Karney (2003) is used. The other two
  methods are faster but less precise

- pairs:

  logical. If `TRUE` pairwise nearest points are returned (only relevant
  when using at least one SpatVector of lines or polygons

- lines:

  logical. If `TRUE` lines between the nearest points instead of (the
  nearest) points

## See also

[`distance`](https://rspatial.github.io/terra/reference/distance.md),
[`relate`](https://rspatial.github.io/terra/reference/relate.md),
[`adjacent`](https://rspatial.github.io/terra/reference/adjacent.md)

## Value

matrix

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
nearby(v, distance=12000)
#>      from to
#> [1,]    2  4
#> [2,]    6  8
```
