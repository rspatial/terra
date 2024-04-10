\docType{methods}
\name{flowAccu2_weight}

\alias{flowAccu2_weight}
\alias{flowAccu2_weight,SpatRaster,SpatRaster-method}


\title{Flow accumulation or Total contributing Area  (Watershed Analysis)}

\description{
It computes Flow accumulation or Total contributing Area (Waterhed Analysis) or a cumulation value in function of a weight value asspciated to each cell. 
}

\usage{
\S4method{flowAccu2_weight}{SpatRaster,SpatRaster}(p,weight,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{weight}{Weight/Score raster map, a \code{\link{SpatRaster-class}} object,  e.g call area map or Averaged Precipitation Map. }
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing the value for each cell. 
}

\details{
The algoritm is an adaptation of the one proposed by Zhou et al, 2019


Respect to the one of \code{\link{flowAccu2}}, this function  takes into account that each cell can contain a value/weight/score, see \code{weight} argument, e.g. the planar area of each cell (in square kilometres) or a precipitation runoff value, in this case the flow accumulations is the total contributing area or the total contributing runoff W. In \code{\link{flowAccu2}} ,instead, each cell contains a value equal to 1 and  the flow accumulation corresponds only to the number of upslope cells. 

Further ArGIS implementation of this index: \url{https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-accumulation.html}
}

\seealso{\code{\link{flowAccu2}},\code{\link{NIDP2}}}
\author{

Emanuele Cordano
}

\examples{


##source("~/local/rpackages/jrc/terra/R/generics_watershed2_extension.R")


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
##
weight <- elev1*0+10

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

flow_acc1 <- flowAccu2(flowdir1)
flow_acc2 <- flowAccu2(flowdir2)

t(array(flow_acc1[],rev(dim(flow_acc1)[1:2])))
t(array(flow_acc2[],rev(dim(flow_acc2)[1:2])))

plot(flow_acc1)
plot(flow_acc2)




flow_acc1w <- flowAccu2_weight(flowdir1,weight)
flow_acc2w <- flowAccu2_weight(flowdir2,weight)

t(array(flow_acc1w[],rev(dim(flow_acc1w)[1:2])))
t(array(flow_acc2w[],rev(dim(flow_acc2w)[1:2])))

plot(flow_acc1w)
plot(flow_acc2w)


## Application wth example DEM

elev <- rast(system.file('ex/elev.tif',package="terra"))
flowdir <- terrain(elev,"flowdir")

weight <- cellSize(elev,unit="km")
flowacc_weight <- flowAccu2_weight(flowdir,weight)
flowacc  <- flowAccu2(flowdir)




}




\references{

Zhou, G., Wei, H. & Fu, S. A fast and simple algorithm for calculating flow accumulation matrices from raster digital elevation. Front. Earth Sci. 13, 317–326 (2019). https://doi.org/10.1007/s11707-018-0725-9
\url{https://link.springer.com/article/10.1007/s11707-018-0725-9}
Further references: \url{https://ica-abs.copernicus.org/articles/1/434/2019/}

}



\keyword{spatial}
