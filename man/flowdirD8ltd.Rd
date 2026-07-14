\docType{methods}
\name{flowdirD8ltd}

\alias{flowdirD8ltd}
\alias{flowdirD8ltd,SpatRaster-method}
\alias{flowdirD8lad}
\alias{flowdirD8lad,SpatRaster-method}

\title{Path-Based Nondisperive Flow Direction}

\description{
Computes non-disperive flow/draiage direction according to path-based methods over grid-based digital elevation models. 
}

\usage{
\S4method{flowdirD8ltd}{SpatRaster}(x,lambda=0.5,deviation_type=c("ltd","lad"),max_iters=10^6,filename="", ...) 
\S4method{flowdirD8lad}{SpatRaster}(x,lambda=0.5,deviation_type=c("lad","ltd"),max_iters=10^6,filename="",...) 
}



\arguments{
  \item{x}{SpatRaster with elevation/terrain model , see \code{\link{terrain}}. }
  \item{lambda}{Parameter between 0 and 1.}  
  \item{deviation_type}{Character string. Default is the first element of \code{=c("ltd","lad")}. If \code{"ltd"} (default) flow direction dispersion is deteceted with LTD criterion, if \code{"lad"} flow direction dispersion is deteceted with LTD criterion. See Orlandini et al,2003 for further details.  }  
  \item{max_iters}{maximum iterations for drainage path starting points detection}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster 
}

\details{
The algorithm is an adaptation of the one proposed by Li et al, 2021 and Orlandini et al,2003
This function is experimental and under development: results are to be verified. 
}

\seealso{\code{\link{terrain}},\code{\link{watershed}},\code{\link{flowAccumulation}}}

\author{
Zhenya Li at al, Stefano Orlandini at al. , Emanuele Cordano (R/C/C++ code implemetation)
}



\references{

Orlandini, S., G. Moretti, M. Franchini, B. Aldighieri, and B. Testa (2003), Path-based methods for the determination of nondispersive drainage directions in grid-based digital elevation models, Water Resour. Res., 39, 1144, doi:10.1029/2002WR001639, 6. \url{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002WR001639}

Li, Z., Shi, P., Yang, T., Wang, C., Yong, B., & Song, Y. (2022). An improved D8-LTD for the extraction of total contributing area (TCA) by adopting the strategies of path independency and local dispersion. Water Resources Research, 58, e2021WR030948. https://doi.org/10.1029/2021WR030948 \url{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021WR030948}


Useful presentation: \url{http://www.idrologia.unimore.it/orlandini/web-archive/seminars/nyc-2008-2.pdf}

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
flowdir1l1<- flowdirD8ltd(elev1,lambda=1)
flowdir2l1<- flowdirD8ltd(elev2,lambda=1)

t(array(flowdir1l1[],rev(dim(flowdir1l1)[1:2])))
t(array(flowdir2l1[],rev(dim(flowdir2l1)[1:2])))


flowdir1ldf<- flowdirD8ltd(elev1,lambda=0.5)
flowdir2ldf<- flowdirD8ltd(elev2,lambda=0.5)

t(array(flowdir1ldf[],rev(dim(flowdir1l1)[1:2])))
t(array(flowdir2ldf[],rev(dim(flowdir2l1)[1:2])))




flowdir1l0<- flowdirD8ltd(elev1,lambda=0)
flowdir2l0<- flowdirD8ltd(elev2,lambda=0)

t(array(flowdir1l0[],rev(dim(flowdir1l1)[1:2])))
t(array(flowdir2l0[],rev(dim(flowdir2l1)[1:2])))


## Flow direction raster
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


## Cone Geometry using calculus R package 


library(calculus)
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "(x^2+y^2)^(1/2)" ## cone ##ok
###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x \%in\% range(vars$x))] <- NA 
vars$elev[which(vars$y \%in\% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0)
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5)
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1)

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")

## V-Shape Watershed  Geometry using calculus R package 

library(calculus)
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "x*tanh(100*x)+0.2*y" ## v-shape basin
###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x \%in\% range(vars$x))] <- NA 
vars$elev[which(vars$y \%in\% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirladl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0)
rr$flow_dirladlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0.5)
rr$flow_dirladl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=1)
rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="ltd")
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="ltd")
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="ltd")

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl0, unit="flowdir",col="blue")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladlm, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl1, unit="flowdir",col="black")

plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")



## Planar Hillslope  Geometry using calculus R package 

library(calculus)
library(stringr)

dx <- 1 #2.5
dy <- 1 #2.5
xmin <- -10
ymin <- -10

xmax <- 10
ymax <- 10


x <- seq(from=xmin,to=xmax,by=dx)
y <- seq(from=ymin,to=ymax,by=dy)


####
vars <- list()
vars$x <- rep(x,times=length(y))
vars$y <- rep(y,each=length(x))
vars <- as.data.frame(vars)
### 
elev_f <- "x+0.5*y" ## planar hillslope basin , try also with elev_f <- "x+0.505*y"

###

vars$elev <- evaluate(elev_f,var=vars[,c("x","y")])
vars$elev[which(vars$x \%in\% range(vars$x))] <- NA 
vars$elev[which(vars$y \%in\% range(vars$y))] <- NA 

## Flow Angle 

vgrad <- matrix(gradient(elev_f,var=vars[,c("x","y")],accuracy=8),ncol=2)
vars$flow_angle <- atan2(y=-vgrad[,2],x=-vgrad[,1])+vars$elev*0


####


rr <- rast(vars)
rr$flow_dir <- terrain(rr$elev,"flowdir")
rr$flow_dirladl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0) 

rr$flow_dirladlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="lad") 
##or flowdirD8lad(rr$elev,lambda=0.5) 

rr$flow_dirladl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="lad")  
##or flowdirD8lad(rr$elev,lambda=1) 

rr$flow_dirltdl0 <- flowdirD8ltd(rr$elev,lambda=0,deviation_type="ltd")
rr$flow_dirltdlm <- flowdirD8ltd(rr$elev,lambda=0.5,deviation_type="ltd")
rr$flow_dirltdl1 <- flowdirD8ltd(rr$elev,lambda=1,deviation_type="ltd")

plot(rr$elev)
arrows_on_rast(rr$flow_dir, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl0, unit="flowdir",col="blue")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladlm, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirladl1, unit="flowdir",col="black")

plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl0, unit="flowdir",col="blue")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdlm, unit="flowdir",col="black")
plot(rr$elev)
arrows_on_rast(rr$flow_dirltdl1, unit="flowdir",col="black")





### Further examples

elev <- rast(system.file('ex/elev.tif',package="terra"))

flowdirlad <- flowdirD8lad(elev,lambda=0.5)
flowdirlad_copy <- flowdirD8ltd(elev,lambda=0.5,deviation_type="lad")

if (!all(values(flowdirlad_copy==flowdirlad),na.rm=TRUE)) {
  stop("Something in flowdirD8lad  went wrong!") 
}

flowdirladl0 <- flowdirD8lad(elev,lambda=0)
flowdirladl0_copy <- flowdirD8ltd(elev,lambda=0,deviation_type="lad")

if (!all(values(flowdirladl0_copy==flowdirladl0),na.rm=TRUE)) {
  stop("Something in flowdirD8lad  went wrong!") 
}

flowdir <- terrain(elev,"flowdir")
flowdir[flowdirladl0==0] <- 0 
flowdir==flowdirladl0




}

\keyword{spatial}
