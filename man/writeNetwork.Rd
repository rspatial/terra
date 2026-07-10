\name{writeNetwork}

\docType{methods}

\alias{writeNetwork}
\alias{writeNetwork,SpatNetwork,character-method}

\title{Write a SpatNetwork to disk}

\description{
Write a \code{\link{SpatNetwork}} to disk as a GDAL Geographic Network Model (GNM) dataset. Use \code{\link{netw}(filename)} to read it back in.

A GNM dataset stores the network topology, the node and edge geometries, the edge weights, the directedness flag and the network's spatial reference in a self-contained directory (\code{filetype="GNMFile"}) or single OGR-compatible database (\code{filetype="GNMDatabase"}). It is the only widely-supported on-disk format that round-trips a topological network with both geometry and attributes.
}

\usage{
\S4method{writeNetwork}{SpatNetwork,character}(x, filename, filetype="GNMFile", overwrite=FALSE, options=NULL, ...)
}

\arguments{
\item{x}{\code{SpatNetwork}.}
\item{filename}{character. Path to write to. For \code{filetype="GNMFile"} this is the path of the network directory; the directory must not exist (or \code{overwrite=TRUE} must be passed). For \code{filetype="GNMDatabase"} this is the connection string of the underlying database.}
\item{filetype}{character. One of \code{"GNMFile"} (the default; a directory of small OGR layer files) or \code{"GNMDatabase"} (a single OGR-compatible database).}
\item{overwrite}{logical. If \code{TRUE} and \code{filename} already exists, it is removed before writing.}
\item{options}{character. Optional driver-specific creation options, formatted as \code{"NAME=VALUE"} (see the GDAL GNM driver documentation). Defaults to no options.}
\item{...}{additional arguments (currently ignored).}
}

\details{
A \code{"GNMFile"} dataset is a directory containing two "class" layers (\code{nodes} and \code{edges}) plus the GNM system layers (\code{_gnm_meta}, \code{_gnm_graph}, \code{_gnm_features}, \code{_gnm_srs.prj}). The default backend is ESRI Shapefile, so a freshly-written network is just a directory of \code{.shp}/\code{.dbf}/\code{.shx}/\code{.prj}/\code{.dbf} files plus the system files.

What round-trips through GNM and what doesn't:
\itemize{
\item Preserved exactly: node coordinates, edge geometries, the topology incidence list, the edge length cache, the edge weight cache, directedness, and the network's CRS.
\item Lost: \code{source_id} attribution and any user-attached columns on \code{net_edges()} or \code{net_nodes()} that go beyond the GNM-defined fields. (GNM has no concept of user attribute schemas on its class layers.)
}

GNM support requires GDAL \eqn{\geq} 2.4 built with the GNM component enabled. If your GDAL was built without GNM, this function fails with a message and \code{\link{netw}(filename)} cannot read GNM datasets back in.
}

\value{
Invisibly \code{TRUE} on success. On failure terra raises an error.
}

\seealso{
\code{\link{netw}}, \code{\link{shortestPath}}
}

\examples{
\dontrun{
v <- vect(rbind(
  "LINESTRING(0 0, 10 10)",
  "LINESTRING(0 10, 10 0)"
), crs = "EPSG:32633")
n <- netw(v)

dst <- file.path(tempdir(), "mynet")
writeNetwork(n, dst, overwrite = TRUE)

# Read back through netw():
n2 <- netw(dst)
}
}


\keyword{methods}
\keyword{spatial}
