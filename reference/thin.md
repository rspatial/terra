# Subset geometries by minimum distance

`thin` return a subset of geometries such that all remaining geometries
are at least `d` apart from each other. This uses greedy spatial
thinning: the first geometry is always kept, and subsequent geometries
are only kept if they are at least `d` away from all previously kept
geometries.

`thinNodes` reduces the number of nodes to represent line or polygon
geometries

## Usage

``` r
# S4 method for class 'SpatVector'
thin(x, d, unit="m")

# S4 method for class 'SpatVector'
thinNodes(x, d, unit="m", makeValid=TRUE)
```

## Arguments

- x:

  SpatVector

- d:

  positive numeric. The minimum distance between geometries or nodes

- unit:

  character. `"m"` (meter, the default) or `"km"` (kilometer)

- makeValid:

  logical. If `TRUE` an attempt is made to assure that thinned polygon
  geometries are valid (e.g., not self intersecting)

## Value

SpatVector

## See also

[`densify`](https://rspatial.github.io/terra/reference/densify.md),
[`distance`](https://rspatial.github.io/terra/reference/distance.md),
[`nearby`](https://rspatial.github.io/terra/reference/nearby.md),
[`nearest`](https://rspatial.github.io/terra/reference/nearby.md)

## Examples

``` r
p <- vect(rbind(c(0,0), c(0.5,0.5), c(0.6, 0.6), c(5,5), c(5.5, 5.5), c(50, 50)), crs="+proj=longlat +datum=WGS84")
values(p) <- data.frame(id=1:6)
nrow(p)
#> [1] 6

# points that are at least 200 km apart
p2 <- thin(p, 200000)
nrow(p2)
#> [1] 3

# same, using km
p3 <- thin(p, 200, unit="km")
nrow(p3)
#> [1] 3
```
