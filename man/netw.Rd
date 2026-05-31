\name{netw}

\docType{methods}

\alias{netw}
\alias{netw,missing-method}
\alias{netw,SpatVector-method}
\alias{netw,igraph-method}
\alias{netw,character-method}
\alias{netw,SpatNetwork-method}
\alias{net_nodes}
\alias{net_nodes,SpatNetwork-method}
\alias{net_edges}
\alias{net_edges,SpatNetwork-method}
\alias{net_nnodes}
\alias{net_nnodes,SpatNetwork-method}
\alias{net_nedges}
\alias{net_nedges,SpatNetwork-method}
\alias{net_directed}
\alias{net_directed,SpatNetwork-method}
\alias{net_weights}
\alias{net_weights,SpatNetwork-method}
\alias{net_weights<-}
\alias{net_weights<-,SpatNetwork-method}
\alias{plot,SpatNetwork,missing-method}
\alias{show,SpatNetwork-method}
\alias{crs,SpatNetwork-method}
\alias{ext,SpatNetwork-method}
\alias{coerce,SpatVector,SpatNetwork-method}
\alias{coerce,igraph,SpatNetwork-method}
\alias{coerce,SpatNetwork,igraph-method}

\title{Build a SpatNetwork}

\description{
\code{netw} is the user-facing constructor for a \code{SpatNetwork}, in the same spirit as \code{\link{vect}} for \code{SpatVector} or \code{\link{rast}} for \code{SpatRaster}. It dispatches on the class of its first argument:
\itemize{
  \item \code{netw()} -- an empty network.
  \item \code{netw(<SpatVector of lines>, ...)} -- builds a planar topological network by noding the input lines (every place where two lines cross becomes a node, the line pieces between consecutive nodes become edges). The full coordinate sequence between nodes is preserved on each edge so that the original sinuous geometry survives.
  \item \code{netw(<igraph>, ...)} -- builds a network from an \code{igraph} object. Vertex attributes give the node coordinates; an optional \code{weight} edge attribute supplies edge weights. Edges become straight 2-point lines between their end nodes (igraph carries no edge geometry).
  \item \code{netw(<filename>, ...)} -- reads a network from disk. Currently only the GDAL Geographic Network Model (GNM) format is supported. Use \code{\link{writeNetwork}} for the write side.
  \item \code{netw(<SpatNetwork>, ...)} -- identity (returns the input unchanged).
}

The same coercions are also exposed via \code{methods::as}: \code{as(v, "SpatNetwork")}, \code{as(g, "SpatNetwork")}, \code{as(net, "igraph")}.

The accessors \code{net_nodes}, \code{net_edges}, \code{net_nnodes}, \code{net_nedges}, \code{net_directed} and \code{net_weights} return parts of the network as ordinary \code{SpatVector}s and primitives. The \code{net_*} prefix avoids collisions with the \code{igraph} and \code{network} packages, which export plain \code{nodes}, \code{edges}, \code{is.directed} and \code{weights}.
}

\usage{
\S4method{netw}{missing}(x, ...)

\S4method{netw}{SpatVector}(x, snap=0, merge=TRUE, directed=FALSE, weights=TRUE, ...)

\S4method{netw}{igraph}(x, x_attr="x", y_attr="y", weight_attr="weight", crs=NULL, ...)

\S4method{netw}{character}(x, ...)

\S4method{netw}{SpatNetwork}(x, ...)

\S4method{net_nodes}{SpatNetwork}(x, ...)
\S4method{net_edges}{SpatNetwork}(x, ...)
\S4method{net_nnodes}{SpatNetwork}(x, ...)
\S4method{net_nedges}{SpatNetwork}(x, ...)
\S4method{net_directed}{SpatNetwork}(x, ...)
\S4method{net_weights}{SpatNetwork}(x, ...)
}

\arguments{
\item{x}{the input. For \code{netw}: a \code{SpatVector} of lines, an \code{igraph}, a filename pointing to a GNM dataset, an existing \code{SpatNetwork}, or missing (for an empty network). For the accessor methods: a \code{SpatNetwork}.}
\item{snap}{numeric. Snapping tolerance, in map units, for \code{netw,SpatVector-method}. If \code{> 0}, line endpoints (and other vertices) within \code{snap} of each other are coalesced before noding. Use this to repair small digitizing gaps in road centerlines. Default \code{0}.}
\item{merge}{logical. If \code{TRUE} (the default), runs of degree-2 shared endpoints are stitched back into single edges, so that only true junctions and dangling endpoints become nodes. Set \code{FALSE} to keep every endpoint of every input feature as a node (useful when the input has been pre-segmented and you want to preserve those breaks).}
\item{directed}{logical. If \code{TRUE}, the network is directed and the orientation of each edge (\code{from_node} \eqn{\to} \code{to_node}) is the legal travel direction. Use this for water flow or one-way streets. Default \code{FALSE}.}
\item{weights}{One of: \code{TRUE} (the default) to weight each edge by its geometric length (in meters when the CRS is lon/lat, otherwise in CRS units scaled to meters when the unit is known); \code{FALSE} for an unweighted network; a numeric vector of length \code{nedges} for custom edge weights; or a single character string giving the name of a numeric column of \code{x} -- in which case each output edge gets the column value of its source feature (with edges that could not be attributed falling back to their geometric length, with a warning).}
\item{x_attr, y_attr}{For \code{netw,igraph-method}: names of the vertex attributes giving node coordinates. Defaults match what \code{as(net, "igraph")} writes.}
\item{weight_attr}{For \code{netw,igraph-method}: name of the edge attribute to use as weights. Defaults to \code{"weight"} (the igraph convention).}
\item{crs}{character. For \code{netw,igraph-method}: optional CRS to assign to the \code{SpatNetwork}. If \code{NULL}, the \code{crs} graph attribute is used when present.}
\item{...}{additional arguments (currently ignored).}
}

\details{
The construction from a \code{SpatVector} uses GEOS to "node" the input: all input lines are merged and split at every interior crossing, and exactly-coincident pieces are dissolved. With \code{merge=TRUE}, runs of edges connected through degree-2 nodes are then re-fused into a single edge -- so a long road that was digitized as many short segments shows up as one edge per stretch between real junctions.

Per-edge attribute forwarding: each output edge is associated with the first input feature whose geometry covers the edge midpoint (within \code{snap}, or a small fraction of the bounding box if \code{snap == 0}). The 1-based row index of that input feature is stored as the \code{source_id} column on \code{net_edges()}, and the corresponding row of the input \code{data.frame} is copied into the edge attribute table. When \code{merge=TRUE}, an edge can span multiple source features that happened to be collinear and connected; only one source's attributes are forwarded in that case.

\code{net_nodes()} returns a \code{SpatVector} of points, one per node, with a \code{degree} column (number of incident edges). For directed networks, \code{in_degree} and \code{out_degree} columns are added.

\code{net_edges()} returns a \code{SpatVector} of lines, one per edge, with columns \code{from_node}, \code{to_node}, \code{source_id}, \code{length} (cached geometric length, in meters for lon/lat) and -- when the network is weighted -- a \code{weight} column. Any forwarded source attributes follow.

\code{net_weights(net)} returns the per-edge weights (or \code{NULL} for unweighted networks). \code{net_weights(net) <- v} replaces the weights with a numeric vector of length \code{net_nedges(net)}; setting it to \code{NULL} marks the network as unweighted.

For the \code{igraph} entry point: forward (\code{SpatNetwork -> igraph}) preserves topology, directedness, edge weights, node coordinates (as vertex attributes \code{x} and \code{y}), additional edge columns from \code{net_edges()}, and the CRS (WKT2, on the graph attribute \code{crs}). Reverse (\code{igraph -> SpatNetwork}) is intentionally lossy: each output edge becomes a straight 2-point line between its end nodes since igraph carries no edge geometry. The \code{igraph} package is a Suggested dependency of terra; install it before calling these methods.
}

\value{
A \code{SpatNetwork} for \code{netw}; a \code{SpatVector} for \code{net_nodes}/\code{net_edges}; an integer for \code{net_nnodes}/\code{net_nedges}; logical for \code{net_directed}; numeric (or \code{NULL}) for \code{net_weights}.
}

\seealso{\code{\link{shortestPath}}, \code{\link{writeNetwork}}}

\examples{
# Two crossing roads
a <- vect("LINESTRING(0 0, 10 10)", crs="local")
b <- vect("LINESTRING(0 10, 10 0)", crs="local")
roads <- rbind(a, b)
roads$name <- c("A", "B")

net <- netw(roads)
net
net_nnodes(net)   # 5: 4 ends + 1 crossing
net_nedges(net)   # 4: each road is split in two
net_weights(net)  # length-based, all sqrt(50)

# Directed network (e.g. a stream network: from_node -> to_node is downstream)
streams <- netw(roads, directed=TRUE)
net_directed(streams)
net_nodes(streams)  # in_degree, out_degree columns appear

# Custom weights
net_weights(net) <- net_weights(net) * 1.5

# Coercion via setAs
net2 <- as(roads, "SpatNetwork")

\dontrun{
# igraph round-trip (requires igraph)
g    <- as(net, "igraph")
back <- netw(g)        # or as(g, "SpatNetwork")
}

# Plot
plot(net)
}


\keyword{methods}
\keyword{spatial}
