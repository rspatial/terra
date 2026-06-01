# Write a SpatNetwork to disk

Write a
[`SpatNetwork`](https://rspatial.github.io/terra/reference/netw.md) to
disk as a GDAL Geographic Network Model (GNM) dataset. Use
[`netw`](https://rspatial.github.io/terra/reference/netw.md)`(filename)`
to read it back in.

A GNM dataset stores the network topology, the node and edge geometries,
the edge weights, the directedness flag and the network's spatial
reference in a self-contained directory (`filetype="GNMFile"`) or single
OGR-compatible database (`filetype="GNMDatabase"`). It is the only
widely-supported on-disk format that round-trips a topological network
with both geometry and attributes.

## Usage

``` r
# S4 method for class 'SpatNetwork,character'
writeNetwork(x, filename, filetype="GNMFile", overwrite=FALSE, options=NULL, ...)
```

## Arguments

- x:

  `SpatNetwork`.

- filename:

  character. Path to write to. For `filetype="GNMFile"` this is the path
  of the network directory; the directory must not exist (or
  `overwrite=TRUE` must be passed). For `filetype="GNMDatabase"` this is
  the connection string of the underlying database.

- filetype:

  character. One of `"GNMFile"` (the default; a directory of small OGR
  layer files) or `"GNMDatabase"` (a single OGR-compatible database).

- overwrite:

  logical. If `TRUE` and `filename` already exists, it is removed before
  writing.

- options:

  character. Optional driver-specific creation options, formatted as
  `"NAME=VALUE"` (see the GDAL GNM driver documentation). Defaults to no
  options.

- ...:

  additional arguments (currently ignored).

## Details

A `"GNMFile"` dataset is a directory containing two "class" layers
(`nodes` and `edges`) plus the GNM system layers (`_gnm_meta`,
`_gnm_graph`, `_gnm_features`, `_gnm_srs.prj`). The default backend is
ESRI Shapefile, so a freshly-written network is just a directory of
`.shp`/`.dbf`/`.shx`/`.prj`/`.dbf` files plus the system files.

What round-trips through GNM and what doesn't:

- Preserved exactly: node coordinates, edge geometries, the topology
  incidence list, the edge length cache, the edge weight cache,
  directedness, and the network's CRS.

- Lost: `source_id` attribution and any user-attached columns on
  [`net_edges()`](https://rspatial.github.io/terra/reference/netw.md) or
  [`net_nodes()`](https://rspatial.github.io/terra/reference/netw.md)
  that go beyond the GNM-defined fields. (GNM has no concept of user
  attribute schemas on its class layers.)

GNM support requires GDAL \\\geq\\ 2.4 built with the GNM component
enabled. If your GDAL was built without GNM, this function fails with a
message and
[`netw`](https://rspatial.github.io/terra/reference/netw.md)`(filename)`
cannot read GNM datasets back in.

## Value

Invisibly `TRUE` on success. On failure terra raises an error.

## See also

[`netw`](https://rspatial.github.io/terra/reference/netw.md),
[`shortestPath`](https://rspatial.github.io/terra/reference/shortestPath.md)

## Examples

``` r
if (FALSE) { # \dontrun{
v <- vect(rbind(
  "LINESTRING(0 0, 10 10)",
  "LINESTRING(0 10, 10 0)"
), crs = "EPSG:32633")
n <- netw(v)

dst <- file.path(tempdir(), "mynet")
writeNetwork(n, dst, overwrite = TRUE)

# Read back through netw():
n2 <- netw(dst)
} # }
```
