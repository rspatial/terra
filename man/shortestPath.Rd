\name{shortestPath}

\docType{methods}

\alias{shortestPath}
\alias{shortestPath,SpatNetwork-method}

\title{Shortest paths on a SpatNetwork}

\description{
Compute shortest paths between pairs of nodes on a \code{SpatNetwork} using Dijkstra's algorithm. The path weight is the sum of \code{\link{net_weights}} along the traversed edges -- by default this is geometric length in meters (geodesic for lon/lat data, planar Euclidean otherwise), but custom weights such as travel time can be set with \code{net_weights<-} or via the \code{weights} argument of \code{\link{netw}}.

The function honours \code{net_directed(net)}: in a directed network only edges traversed in their stored \code{from_node \\eqn{\\to} to_node} orientation are allowed, which is what you want for one-way streets or stream networks.

The result is a \code{SpatVector} of lines that follow the original (sinuous) edge geometries, so that the returned features can be plotted, projected and written like any other \code{SpatVector}.
}

\usage{
\S4method{shortestPath}{SpatNetwork}(x, from, to, ...)
}

\arguments{
\item{x}{\code{SpatNetwork}}
\item{from}{integer. One or more 1-based node ids (rows of \code{net_nodes(x)}) to start from.}
\item{to}{integer. One or more 1-based node ids to reach. Either \code{from} or \code{to} may be length 1 to be recycled against the other; otherwise the two must have the same length.}
\item{...}{additional arguments (currently ignored)}
}

\details{
The implementation builds a single adjacency list from \code{net_edges(x)} and then runs one Dijkstra single-source shortest-path traversal per unique source node, caching the result so that querying many destinations from the same origin is efficient. Negative or \code{NA} edge weights are treated as zero (Dijkstra cannot handle negative cycles).

When a destination is unreachable from its source the corresponding output feature has an empty geometry and \code{NA} distance. When \code{from == to} the distance is \code{0} and the geometry is empty.
}

\value{
A \code{SpatVector} of lines with one feature per (\code{from}, \code{to}) pair, with columns
\itemize{
\item \code{from}: 1-based source node id
\item \code{to}: 1-based destination node id
\item \code{distance}: total path weight (\code{NA} if unreachable)
}
}

\seealso{\code{\link{netw}}, \code{\link{net_weights}}}

\examples{
# Two crossing roads
a <- vect("LINESTRING(0 0, 10 10)", crs = "local")
b <- vect("LINESTRING(0 10, 10 0)", crs = "local")
roads <- rbind(a, b)
roads$speed <- c(50, 100)   # km/h, used as a "cost" column

net <- netw(roads)
sp <- shortestPath(net, from = 1, to = 4)
sp$distance

# Multiple destinations from one source
shortestPath(net, from = 1, to = c(2, 3, 4, 5))$distance

# Custom weights from a column of the source SpatVector
nslow <- netw(roads, weights = "speed")
shortestPath(nslow, from = 1, to = 4)$distance
}


\keyword{methods}
\keyword{spatial}
