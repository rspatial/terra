\docType{methods}
\name{watershed2}

\alias{watershed2}
\alias{watershed2,SpatRaster,integer-method}


\title{Watershed Extraction}

\description{
It extracts a watershed mask from a flow direction input raster map and  a pour-point cell (watershed outlet).
}

\usage{
\S4method{watershed2}{SpatRaster,integer}(p,pp_offset,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{pp_offset}{pour point offset.}
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 in the extracted watershed mask and value 0 elsewhere. 
}

\details{

This algorithm was created and tested by the authors. 
}


\author{

Ezio Crestaz, Emanuele Cordano, Roman Seliger
}

\examples{

library(sf)



elev <- rast(system.file('ex/elev.tif',package="terra"))
flowdir <- terrain(elev,"flowdir")
flowacc  <- flowAccu2(flowdir)
## pourpoint with maximum accumulation. 
pp <- which.max(flowacc[])


poutr <- watershed2((flowdir),as.integer(pp))

poutr[poutr==0] <- NA
poutv <- st_as_sf(as.polygons(poutr))
plot(st_geometry(poutv))
plot(elev+poutr*0,add=TRUE)



}




\keyword{spatial}
