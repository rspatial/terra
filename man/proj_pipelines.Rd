\name{proj_pipelines}

\alias{proj_pipelines}

\title{Find CRS transformation pipelines}

\description{
Find candidate coordinate transformation pipelines between two coordinate reference systems. This can be used to examine which transformations are available, whether they require grid files, and what accuracy they provide.
}

\usage{
proj_pipelines(from, to, authority="", AOI=NULL, use="NONE", grid_availability="USED", 
	desired_accuracy=-1.0, strict_containment=FALSE, axis_order_authority_compliant=FALSE)
}

\arguments{
\item{from}{character, SpatRaster, SpatVector, or any object from which \code{crs()} can extract a CRS}
\item{to}{same types as for \code{from}}
\item{authority}{character. Constrain output pipelines to those from this authority (e.g. \code{"EPSG"}). Default \code{""} means no constraint}
\item{AOI}{SpatExtent. Area of interest for the transformation in degrees. For an area of interest crossing the anti-meridian, west should be greater than east, e.g. \code{ext(150, 210, 0, 30)}}
\item{use}{character. One of \code{"NONE"}, \code{"BOTH"}, \code{"INTERSECTION"}, \code{"SMALLEST"}, indicating how areas of interest of source and target CRS are used. Ignored when \code{AOI} is specified}
\item{grid_availability}{character. One of \code{"USED"} (grid availability is used for sorting; operations with missing grids are sorted last), \code{"DISCARD"} (discard operations if a required grid is missing), \code{"IGNORED"} (ignore grid availability; results presented as if all grids were available), or \code{"AVAILABLE"} (results presented as if grids known to PROJ were available; typical when networking is enabled)}
\item{desired_accuracy}{numeric. Only return pipelines with at least this accuracy (in metres). Use \code{-1} (the default) to impose no accuracy constraint}
\item{strict_containment}{logical. If \code{FALSE} (the default), permit partial matching of the area of interest. If \code{TRUE}, the area of interest must be strictly contained}
\item{axis_order_authority_compliant}{logical. If \code{FALSE} (the default), always use longitude/easting for the first axis. If \code{TRUE}, follow the axis order given by the CRS definitions}
}

\value{
A \code{data.frame} with "pipelines" and the following columns:
\describe{
  \item{description}{human-readable description}
  \item{definition}{PROJ pipeline definition string}
  \item{has_inverse}{logical; whether the pipeline has an inverse}
  \item{accuracy}{accuracy in metres; \code{-1} indicates unknown (ballpark) accuracy}
  \item{grid_count}{number of grids required by this pipeline}
  \item{instantiable}{logical; whether the pipeline can be used (all required grids are available)}
}
}

\note{
Requires PROJ >= 7.1 (inc" 3). See \code{\link{projNetwork}} to enable online grid downloading.
}

\seealso{
\code{\link{crs}}, \code{\link{project}}, \code{\link{same.crs}}, \code{\link{projNetwork}}
}

\examples{
# Simple case: WGS84 to UTM zone 32Nproj_pipelines("EPSG:4326", "EPSG:32632")

# Many pipelines between NAD83/Conus Albers and WGS84
pp <- proj_pipelines("EPSG:5070", "EPSG:4326")
nrow(pp)

# Only pipelines that work without downloading grids
pp_local <- proj_pipelines("EPSG:5070", "EPSG:4326", grid_availability="DISCARD")
pp_local

# Using a SpatRaster as input
r <- rast(ncols=10, nrows=10, crs="EPSG:4326")
proj_pipelines(r, "EPSG:3857")
}

\keyword{spatial}
