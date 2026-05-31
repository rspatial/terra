# Build a SpatNetwork from lines

`as.network` builds a planar topological network (a `SpatNetwork`) from
a `SpatVector` of lines, e.g. a road network. Every place where two
input lines cross becomes a node, and the line pieces between
consecutive nodes become edges. Each edge keeps the full coordinate
sequence of the source line between its two end nodes, so the original
sinuous geometry is preserved.

The result is a topological graph (with integer `from`/`to` node indices
per edge) that also retains its full geometry. `nodes()` and `edges()`
return the network as ordinary `SpatVector` objects so that all of
terra's plotting, projection and writing infrastructure works on them
directly.

## Usage

``` r
# S4 method for class 'SpatVector'
as.network(x, snap=0, merge=TRUE)

# S4 method for class 'SpatNetwork'
nodes(x, ...)
# S4 method for class 'SpatNetwork'
edges(x, ...)
# S4 method for class 'SpatNetwork'
nnodes(x, ...)
# S4 method for class 'SpatNetwork'
nedges(x, ...)
```

## Arguments

- x:

  For `as.network`: a `SpatVector` with `geomtype` `"lines"`. For the
  other methods: a `SpatNetwork`.

- snap:

  numeric. Snapping tolerance, in map units. If `> 0`, line endpoints
  (and other vertices) within `snap` of each other are coalesced before
  noding. Use this to repair small digitizing gaps in road centerlines.
  Default `0`.

- merge:

  logical. If `TRUE` (the default), runs of degree-2 shared endpoints
  are stitched back into single edges, so that only true junctions and
  dangling endpoints become nodes. Set `FALSE` to keep every endpoint of
  every input feature as a node (useful when the input has been
  pre-segmented and you want to preserve those breaks).

- ...:

  additional arguments (currently ignored)

## Details

The construction uses GEOS to "node" the input: all input lines are
merged and split at every interior crossing, and exactly-coincident
pieces are dissolved. With `merge=TRUE`, runs of edges connected through
degree-2 nodes are then re-fused into a single edge – so a long road
that was digitized as many short segments shows up as one edge per
stretch between real junctions.

Per-edge attribute forwarding: each output edge is associated with the
first input feature whose geometry covers the edge midpoint (within
`snap`, or a small fraction of the bounding box if `snap == 0`). The
1-based row index of that input feature is stored as the `source_id`
column on `edges()`, and the corresponding row of the input `data.frame`
is copied into the edge attribute table. When `merge=TRUE`, an edge can
span multiple source features that happened to be collinear and
connected; only one source's attributes are forwarded in that case.

`nodes()` returns a `SpatVector` of points, one per node, with a
`degree` column (number of incident edges).

`edges()` returns a `SpatVector` of lines, one per edge, with columns
`from_node`, `to_node`, `source_id`, `length` (planar Euclidean), and
any forwarded source attributes.

## Value

A `SpatNetwork` for `as.network`; a `SpatVector` for `nodes`/`edges`; an
integer for `nnodes`/`nedges`.

## Examples

``` r
# Two crossing roads
a <- vect("LINESTRING(0 0, 10 10)", crs="local")
b <- vect("LINESTRING(0 10, 10 0)", crs="local")
roads <- rbind(a, b)
roads$name <- c("A", "B")

net <- as.network(roads)
net
#> class       : SpatNetwork
#> dimensions  : 5, 4  (nodes, edges)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter)
#> names       :  name
#> type        : <chr>
nnodes(net)   # 5: 4 ends + 1 crossing
#> [1] 5
nedges(net)   # 4: each road is split in two
#> [1] 4

# Inspect topology
edges(net)
#> class       : SpatVector
#> geometry    : lines
#> dimensions  : 4, 5  (geometries, attributes)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter)
#> names       : from_node to_node source_id  length  name
#> type        :     <int>   <int>     <int>   <num> <chr>
#> values      :         1       2         1 7.07107     A
#>                       3       2         2 7.07107     B
#>                       2       4         1 7.07107     A
#>               ...
nodes(net)
#> class       : SpatVector
#> geometry    : points
#> dimensions  : 5, 1  (geometries, attributes)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter)
#> names       : degree
#> type        :  <int>
#> values      :      1
#>                    4
#>                    1
#>               ...

# Plot
plot(net)
```
