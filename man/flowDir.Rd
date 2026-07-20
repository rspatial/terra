\docType{methods}

\name{flowdir}

\alias{flowDir}
\alias{flowDir,SpatRaster-method}

\title{Path-Based Nondisperive Flow Direction}

\description{
Computes non-disperive flow/draiage direction according to path-based methods over grid-based digital elevation models. 
This is an alternative to \code{\link{terrain}(v="flowdir")}.
}

\usage{
\S4method{flowDir}{SpatRaster}(x, lambda=0.5, deviation_type=c("ltd","lad"), max_iters=10^6, filename="", ...) 
}

\arguments{
  \item{x}{SpatRaster with elevation data}
  \item{lambda}{Parameter between 0 and 1}  
<<<<<<< HEAD
  \item{deviation_type}{Character. Default is the first element of \code{=c("ltd","lad")}. If \code{"ltd"} (default) flow direction dispersion is deteceted with Least Transverse Deviation (LTD) criterion, if \code{"lad"} flow direction dispersion is deteceted with Least Angular Deviation (LAD) criterion. See Orlandini et al., 2003 for details}  
=======
  \item{deviation_type}{Character. Default is the first element of c("ltd", "lad"). If \code{"ltd"} (the default), nondispersive flow directions are determined using the least transversal deviation criterion (LTD). If \code{"lad"}, nondispersive flow directions are determined using the least angular deviation (LAD) criterion. See Orlandini et al. (2003) for details.}  
>>>>>>> 4fb7fade7d313e14147261b00d649a0219e4397d
  \item{max_iters}{maximum iterations for drainage path starting points detection}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster 
}

\details{
The algorithm is an adaptation of the one proposed by Li et al, 2021 and Orlandini et al, 2003
This function is experimental and under development: results are to be verified. 
}

\seealso{\code{\link{terrain}}, \code{\link{watershed}}, \code{\link{flowAccumulation}}}


\author{
Emanuele Cordano
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

## Elevation raster
elev1 <- rast(elev1)
elev2 <- rast(elev2)

t(array(elev1[],rev(dim(elev1)[1:2])))
t(array(elev2[],rev(dim(elev2)[1:2])))

plot(elev1)
plot(elev2)

## Flow direction raster
fdir1 <- flowDir(elev1, lambda=1)
fdir2 <- flowDir(elev2,lambda=1)

elev <- rast(system.file('ex/elev.tif',package="terra"))

fdirlad1 <- flowDir(elev, lambda=0.5, deviation_type="lad")
fdirlad2 <- flowDir(elev, lambda=0.5, deviation_type="lad")
}

\keyword{spatial}
