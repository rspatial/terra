\name{as.network}

\docType{methods}

\alias{as.network}
\alias{as.network,SpatVector-method}
\alias{nodes}
\alias{nodes,SpatNetwork-method}
\alias{edges}
\alias{edges,SpatNetwork-method}
\alias{nnodes}
\alias{nnodes,SpatNetwork-method}
\alias{nedges}
\alias{nedges,SpatNetwork-method}
\alias{plot,SpatNetwork,missing-method}
\alias{show,SpatNetwork-method}
\alias{crs,SpatNetwork-method}
\alias{ext,SpatNetwork-method}

\title{Build a SpatNetwork from lines}

\description{
\code{as.network} builds a planar topological network (a \code{SpatNetwork}) from a \code{SpatVector} of lines, e.g. a road network. Every place where two input lines cross becomes a node, and the line pieces between consecutive nodes become edges. Each edge keeps the full coordinate sequence of the source line between its two end nodes, so the original sinuous geometry is preserved.

The result is a topological graph (with integer \code{from}/\code{to} node indices per edge) that also retains its full geometry. \code{nodes()} and \code{edges()} return the network as ordinary \code{SpatVector} objects so that all of terra's plotting, projection and writing infrastructure works on them directly.
}

\usage{
\S4method{as.network}{SpatVector}(x, snap=0, merge=TRUE)

\S4method{nodes}{SpatNetwork}(x, ...)
\S4method{edges}{SpatNetwork}(x, ...)
\S4method{nnodes}{SpatNetwork}(x, ...)
\S4method{nedges}{SpatNetwork}(x, ...)
}

\arguments{
\item{x}{For \code{as.network}: a \code{SpatVector} with \code{geomtype} \code{"lines"}. For the other methods: a \code{SpatNetwork}.}
\item{snap}{numeric. Snapping tolerance, in map units. If \code{> 0}, line endpoints (and other vertices) within \code{snap} of each other are coalesced before noding. Use this to repair small digitizing gaps in road centerlines. Default \code{0}.}
\item{merge}{logical. If \code{TRUE} (the default), runs of degree-2 shared endpoints are stitched back into single edges, so that only true junctions and dangling endpoints become nodes. Set \code{FALSE} to keep every endpoint of every input feature as a node (useful when the input has been pre-segmented and you want to preserve those breaks).}
\item{...}{additional arguments (currently ignored)}
}

\details{
The construction uses GEOS to "node" the input: all input lines are merged and split at every interior crossing, and exactly-coincident pieces are dissolved. With \code{merge=TRUE}, runs of edges connected through degree-2 nodes are then re-fused into a single edge -- so a long road that was digitized as many short segments shows up as one edge per stretch between real junctions.

Per-edge attribute forwarding: each output edge is associated with the first input feature whose geometry covers the edge midpoint (within \code{snap}, or a small fraction of the bounding box if \code{snap == 0}). The 1-based row index of that input feature is stored as the \code{source_id} column on \code{edges()}, and the corresponding row of the input \code{data.frame} is copied into the edge attribute table. When \code{merge=TRUE}, an edge can span multiple source features that happened to be collinear and connected; only one source's attributes are forwarded in that case.

\code{nodes()} returns a \code{SpatVector} of points, one per node, with a \code{degree} column (number of incident edges).

\code{edges()} returns a \code{SpatVector} of lines, one per edge, with columns \code{from_node}, \code{to_node}, \code{source_id}, \code{length} (planar Euclidean), and any forwarded source attributes.
}

\value{
A \code{SpatNetwork} for \code{as.network}; a \code{SpatVector} for \code{nodes}/\code{edges}; an integer for \code{nnodes}/\code{nedges}.
}

\examples{
# Two crossing roads
a <- vect("LINESTRING(0 0, 10 10)", crs="local")
b <- vect("LINESTRING(0 10, 10 0)", crs="local")
roads <- rbind(a, b)
roads$name <- c("A", "B")

net <- as.network(roads)
net
nnodes(net)   # 5: 4 ends + 1 crossing
nedges(net)   # 4: each road is split in two

# Inspect topology
edges(net)
nodes(net)

# Plot
plot(net)
}


\keyword{methods}
\keyword{spatial}
