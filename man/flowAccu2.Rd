\docType{methods}
\name{flowAccu2}

\alias{flowAccu2}
\alias{flowAccu2,SpatRaster-method}


\title{Flow accumulation or Total contributing Area  (Watershed Analysis)}

\description{
It computes Flow accumulation or Total contributing Area in terms of numbers of cells upstrem each cell(Waterhed Analysis). 
}

\usage{
\S4method{flowAccu2}{SpatRaster}(p,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing the value for each cell. 
}

\details{
The algoritm is an adaptation of the one proposed by Zhou at al, 2019. 



ArGIS impletation (same result but different algorithm):  \url{https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-accumulation.html}
}
\seealso{\code{\link{flowAccu2_weight}},\code{\link{terrain}},\code{\link{watershed2}},\code{\link{NIDP2}}}
\author{

Emanuele Cordano (implementation of the algorithm within the package)
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
    elev2[r,c] <- 10+5*(abs(x))-0.001*y ### 5*(x^2+y^2)
  }
} 


## Elevation Raster Maps
elev1 <- rast(elev1)
elev2 <- rast(elev2)

t(array(elev1[],rev(dim(elev1)[1:2])))
t(array(elev2[],rev(dim(elev2)[1:2])))

plot(elev1)
plot(elev2)

## Flow Direction Raster Maps
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


t(array(flowdir1[],rev(dim(flowdir1)[1:2])))
t(array(flowdir2[],rev(dim(flowdir2)[1:2])))

plot(flowdir1)
plot(flowdir2)

## 
flow_acc1 <- flowAccu2((flowdir1))
flow_acc2 <- flowAccu2((flowdir2))

t(array(flow_acc1[],rev(dim(flow_acc1)[1:2])))
t(array(flow_acc2[],rev(dim(flow_acc2)[1:2])))

plot(flow_acc1)
plot(flow_acc2)


## Application to Central Africa
flowdir_f30 <- rast(system.file('ecor/data/flowdir/nigrb_dir_30s.tif',package="terra")) 
acc_f30 <- flowAccu2(flowdir_f30)
plot(acc_f30>100)


}




\references{
Zhou, G., Wei, H. & Fu, S. A fast and simple algorithm for calculating flow accumulation matrices from raster digital elevation. Front. Earth Sci. 13, 317–326 (2019). https://doi.org/10.1007/s11707-018-0725-9
\url{https://link.springer.com/article/10.1007/s11707-018-0725-9}
Further references: \url{https://ica-abs.copernicus.org/articles/1/434/2019/}
}



\keyword{spatial}
