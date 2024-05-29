\docType{methods}
\name{watershed}

\alias{watershed}
\alias{watershed,SpatRaster-method}


\title{Catchment delineation }

\description{
delineate the area covered b a cachment given from a SpatRaster with flow direction a pour-point (watershed outlet).
}

\usage{
\S4method{watershed}{SpatRaster}(x, pourpoint, filename="",...) 
}

\arguments{
  \item{p}{SpatRaster with flow direction. See \code{\link{terrain}}. }
  \item{pourpoint}{pour point location}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster 
}

\author{
Ezio Crestaz, Emanuele Cordano, Roman Seliger
}

\examples{
elev <- rast(system.file('ex/elev_vinschgau.tif', package="terra"))
flowdir <- terrain(elev, "flowdir")
flowacc <- flowAccu2(flowdir)
## pour point at Naturns 
#pp_lonlat <- vect(x=c(11.004167,46.65)),crs=4326) 
#pp <- projec(pp_lonlat,crs=crs(flowacc))
#ipp <- cellFromXY(flowacc,st_coordinates(pp))
#poutr <- watershed((flowdir),as.integer(ipp))

}

\keyword{spatial}
