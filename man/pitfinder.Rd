\docType{methods}
\name{pitfinder}

\alias{pitfinder}
\alias{pitfinder,SpatRaster-method}


\title{Find pits (depressions with no outlet)}

\description{
find pits (depressions with no outlet) using a flow direction raster, to support the delineation of catchments.
}

\usage{
\S4method{pitfinder}{SpatRaster}(x,pits_on_boundary=TRUE,filename="",...) 
}

\arguments{
  \item{x}{SpatRaster with flow-direcion. See \code{\link{flowDir}}}
   \item{pits_on_boundary}{logical if \code{TRUE} (default) pits are considered also on the boundary. If \code{FALSE} pits are considered only within the terrain domain.}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\value{ 
\code{SpatRaster} with positive integers (1, 2, 3, ...) to identify pits and zero elsewhere. 
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

## Application with example DEM

elev <- rast(system.file('ex/elev.tif',package="terra"))
flowdir <- terrain(elev,"flowdir")

pits <- pitfinder(flowdir)
pits2 <- pitfinder(flowdir,pits_on_boundary=FALSE)
plot((pits>0)==(pits2>0))
}

\keyword{spatial}
