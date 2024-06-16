\docType{methods}
\name{pitfinder}

\alias{pitfinder}
\alias{pitfinder,SpatRaster-method}


\title{Pit Finder in a Flow Dir SpatRaster for Watershed Extraction}

\description{
find pits (depressions with no outlet )
}

\usage{
\S4method{pitfinder}{SpatRaster}(x,filename="",...) 
}

\arguments{
  \item{x}{SpatRaster wih flow-direcion. See \code{\link{terrain}}}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 for the pits  and value 0 elsewhere. 
}


\author{

Emanuele Cordano
}

\seealso{\code{\link{terrain}},\code{\link{watershed}},\code{\link{flowAccumulation}},\code{\link{NIDP}}}
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

pits <- pitfinder(flowdir)

## Application wth example DEM

elev <- rast(system.file('ex/elev.tif',package="terra"))
flowdir <- terrain(elev,"flowdir")

pits <- pitfinder(flowdir)





}




\keyword{spatial}
