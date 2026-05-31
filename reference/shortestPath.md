# Shortest paths on a SpatNetwork

Compute shortest paths between pairs of nodes on a `SpatNetwork` using
Dijkstra's algorithm. The path weight is the sum of
[`net_weights`](https://rspatial.github.io/terra/reference/netw.md)
along the traversed edges – by default this is geometric length in
meters (geodesic for lon/lat data, planar Euclidean otherwise), but
custom weights such as travel time can be set with `net_weights<-` or
via the `weights` argument of
[`netw`](https://rspatial.github.io/terra/reference/netw.md).

The function honours `net_directed(net)`: in a directed network only
edges traversed in their stored `from_node \eqn{\to} to_node`
orientation are allowed, which is what you want for one-way streets or
stream networks.

The result is a `SpatVector` of lines that follow the original (sinuous)
edge geometries, so that the returned features can be plotted, projected
and written like any other `SpatVector`.

## Usage

``` r
# S4 method for class 'SpatNetwork'
shortestPath(x, from, to, ...)
```

## Arguments

- x:

  `SpatNetwork`

- from:

  integer. One or more 1-based node ids (rows of `net_nodes(x)`) to
  start from.

- to:

  integer. One or more 1-based node ids to reach. Either `from` or `to`
  may be length 1 to be recycled against the other; otherwise the two
  must have the same length.

- ...:

  additional arguments (currently ignored)

## Details

The implementation builds a single adjacency list from `net_edges(x)`
and then runs one Dijkstra single-source shortest-path traversal per
unique source node, caching the result so that querying many
destinations from the same origin is efficient. Negative or `NA` edge
weights are treated as zero (Dijkstra cannot handle negative cycles).

When a destination is unreachable from its source the corresponding
output feature has an empty geometry and `NA` distance. When
`from == to` the distance is `0` and the geometry is empty.

## Value

A `SpatVector` of lines with one feature per (`from`, `to`) pair, with
columns

- `from`: 1-based source node id

- `to`: 1-based destination node id

- `distance`: total path weight (`NA` if unreachable)

## See also

[`netw`](https://rspatial.github.io/terra/reference/netw.md),
[`net_weights`](https://rspatial.github.io/terra/reference/netw.md)

## Examples

``` r
# Two crossing roads
a <- vect("LINESTRING(0 0, 10 10)", crs = "local")
b <- vect("LINESTRING(0 10, 10 0)", crs = "local")
roads <- rbind(a, b)
roads$speed <- c(50, 100)   # km/h, used as a "cost" column

net <- netw(roads)
sp <- shortestPath(net, from = 1, to = 4)
sp$distance
#> [1] 14.14214

# Multiple destinations from one source
shortestPath(net, from = 1, to = c(2, 3, 4, 5))$distance
#> [1]  7.071068 14.142136 14.142136 14.142136

# Custom weights from a column of the source SpatVector
nslow <- netw(roads, weights = "speed")
shortestPath(nslow, from = 1, to = 4)$distance
#> [1] 100
```
