\docType{methods}
\name{flowAccumulation}

\alias{flowAccumulation}
\alias{flowAccumulation,SpatRaster-method}


\title{Flow accumulation}

\description{
Computes flow accumulation or the total contributing area in terms of numbers of cells upstream of each cell. 
}

\usage{
\S4method{flowAccumulation}{SpatRaster}(x, weight=NULL, filename="", ...) 
}

\arguments{
  \item{x}{SpatRaster with flow direction, see \code{\link{terrain}}. }
  \item{weight}{SpatRaster with weight/score daa. For example, cell area or precipitation}  
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster
}

\details{
The algorithm is an adaptation of the one proposed by Zhou at al, 2019. 
}

\seealso{\code{\link{terrain}},\code{\link{watershed}},\code{\link{NIDP}}}

\author{
Emanuele Cordano
}


\references{
Zhou, G., Wei, H. & Fu, S. A fast and simple algorithm for calculating flow accumulation matrices from raster digital elevation. Front. Earth Sci. 13, 317–326 (2019). doi:10.1007/s11707-018-0725-9. Also see: \url{https://ica-abs.copernicus.org/articles/1/434/2019/}
}


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
    elev2[r,c] <- 10+5*(abs(x))-0.001*y 
  }
} 


## Elevation Raster Maps
elev1 <- rast(elev1)
elev2 <- rast(elev2)

t(array(elev1[],rev(dim(elev1)[1:2])))
t(array(elev2[],rev(dim(elev2)[1:2])))

plot(elev1)
plot(elev2)

## Flow direction raster
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")
u <- terrain(flowdir2,"flowacc")

t(array(flowdir1[],rev(dim(flowdir1)[1:2])))
t(array(flowdir2[],rev(dim(flowdir2)[1:2])))

plot(flowdir1)
plot(flowdir2)

## 
flow_acc1 <- flowAccumulation((flowdir1))
flow_acc2 <- flowAccumulation((flowdir2))
## Alternative option 
## here the argument of "terrain" functions is  flow directions instead of elevation.

flow_acc2a <- terrain(flowdir2,"flowacc") 
t(array(flow_acc1[],rev(dim(flow_acc1)[1:2])))
t(array(flow_acc2[],rev(dim(flow_acc2)[1:2])))

plot(flow_acc1)
plot(flow_acc2)

weight <- elev1*0+10


flow_acc1w <- flowAccumulation(flowdir1,weight)
flow_acc2w <- flowAccumulation(flowdir2,weight)

t(array(flow_acc1w[],rev(dim(flow_acc1w)[1:2])))
t(array(flow_acc2w[],rev(dim(flow_acc2w)[1:2])))

plot(flow_acc1w)
plot(flow_acc2w)


## Application wth example elevation data

elev <- rast(system.file('ex/elev.tif',package="terra"))
flowdir <- terrain(elev,"flowdir")

weight <- cellSize(elev,unit="km")
flowacc_weight <- flowAccumulation(flowdir,weight)
flowacc  <- flowAccumulation(flowdir)


}

\keyword{spatial}
