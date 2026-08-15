\name{has.geoloc}

\docType{methods}

\alias{has.geoloc}
\alias{has.geoloc,SpatRaster-method}
\alias{geoloc}
\alias{geoloc,SpatRaster-method}


\title{Geolocation arrays and GCPs}

\description{
Detect whether a SpatRaster source has GDAL geolocation arrays (e.g. satellite swath / curvilinear NetCDF) and/or Ground Control Points (GCPs), and report the associated spatial reference and array paths.

When a file has no conventional geotransform but does have a \code{GEOLOCATION} metadata domain (or GCPs), \code{\link{rast}} stores that information on the source and, if the dataset itself has no CRS, applies the geolocation/GCP CRS automatically so that \code{\link{project}} / warping can use the arrays without a manual \code{\link{set.crs}} call.
}

\usage{
\S4method{has.geoloc}{SpatRaster}(x)

\S4method{geoloc}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
}

\value{
\code{has.geoloc}: logical, one value per raster data source.

\code{geoloc}: \code{data.frame} with one row per source and columns
\code{geolocation} (logical), \code{gcps} (logical), \code{srs} (character),
\code{x} / \code{y} (paths to the X/Y geolocation datasets when present).
}

\seealso{
 \code{\link{project}, \link{rectify}, \link{is.rotated}, \link{crs}}
}

\examples{
r <- rast(nrows=10, ncols=10)
has.geoloc(r)
geoloc(r)
}

\keyword{methods}
\keyword{spatial}
