\docType{methods}
\name{pitfinder2}

\alias{pitfinder2}
\alias{pitfinder2,SpatRaster-method}


\title{Pit Finder in a Flow Dir Raster map for Waterhed Extraction}

\description{
It extractions a watershed mask from a flow direction input raster point a pour-point cell.
}

\usage{
\S4method{pitfinder2}{SpatRaster}(p,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 for the pits  and value 0 elsewhere. 
}


\author{

Emanuele Cordano
}

\seealso{\code{\link{terrain}},\code{\link{watershed2}},\code{\link{flowAccu2}},\code{\link{flowAccu2_weight}},\code{\link{NIDP2}}}
\examples{

## Creation of a Digital Elevation Model 

elev <- array(NA,c(9,9))
dx <- 1
dy <- 1 
for (r in 1:nrow(elev)) {
  x <- (r-5)*dx
  for (c in 1:ncol(elev)) {
    
    y <- (c-5)*dy
    elev[r,c] <- 10+5*(x^2+y^2)
    }
  } 
  
elev <- cbind(elev,elev,elev,elev) 
elev <- rbind(elev,elev,elev,elev) 
elev <- rast(elev)

## Flow Directions

flowdir<- terrain(elev,v="flowdir")
t(array(flowdir[],rev(dim(flowdir)[1:2])))

## Pit Detect

pits <- pitfinder2(flowdir)


}




\keyword{spatial}
