\docType{methods}
\name{NIDP}

\alias{NIDP}
\alias{NIDP,SpatRaster-method}


\title{Number of immediate adjacent cells flowing into each cell}

\description{
Compute the number of immediate adjacent cells flowing into each cell
}

\usage{
\S4method{NIDP}{SpatRaster}(x, filename="",...) 
}

\arguments{
  \item{x}{SpatRaster with flow-direction. see \code{\link{terrain}}}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster
}

\details{
NDIP is computed first to compute flow-accumulation with the algorithm by Zhou at al, 2019. 
}




\references{
Zhou, G., Wei, H. & Fu, S. A fast and simple algorithm for calculating flow accumulation matrices from raster digital elevation. Front. Earth Sci. 13, 317–326 (2019). https://doi.org/10.1007/s11707-018-0725-9
\url{https://link.springer.com/article/10.1007/s11707-018-0725-9}
}



\seealso{\code{\link{flowAccumulation}}}

\author{Emanuele Cordano }

\examples{

elev1 <- array(NA,c(9,9))
elev2 <- elev1
dx <- 1
dy <- 1 
for (r in 1:nrow(elev1)) {
  y <- (r-5)*dx
  for (c in 1:ncol(elev1)) {
    
    x <- (c-5)*dy
    elev1[r,c] <- 5*(x^2+y^2)
    elev2[r,c] <- 10+5*(abs(x))-0.001*y ### 5*(x^2+y^2)
  }
} 


## Elevation Raster 
elev1 <- rast(elev1)
elev2 <- rast(elev2)

t(array(elev1[],rev(dim(elev1)[1:2])))
t(array(elev2[],rev(dim(elev2)[1:2])))

plot(elev1)
plot(elev2)

## Flow Direction Raster
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


t(array(flowdir1[],rev(dim(flowdir1)[1:2])))
t(array(flowdir2[],rev(dim(flowdir2)[1:2])))

plot(flowdir1)
plot(flowdir2)

## 
nidp1 <- NIDP((flowdir1))
nidp2 <- NIDP((flowdir2))

t(array(nidp1[],rev(dim(nidp1)[1:2])))
t(array(nidp2[],rev(dim(nidp2)[1:2])))

plot(nidp1)
plot(nidp2)

}



\keyword{spatial}
