\name{overlay}

\docType{methods}

\alias{overlay}
\alias{overlay,SpatRaster,SpatRaster-method}

\title{Overlay SpatRaster objects}

\description{
Create a new SpatRaster object, by combining two or more SpatRaster objects. 

You should supply a function \code{fun} to determine how the SpatRasters are combined. The number of arguments in the function must match the number of SpatRaster objects (or take any number). For example, if you combine two SpatRaster objects you could use multiply: \code{fun=function(x,y){return(x*y)}} percentage: \code{fun=function(x,y){return(100 * x / y)}}. If you combine three layers you could use \code{fun=function(x,y,z){return((x + y) * z)}}

Note that the function work for vectors (not only for single numbers). That is, it must return the same number of elements as its input vectors. 

The function must take each layer as an argument. 

Use \code{\link{app}} for functions such as \code{sum}, that take \code{n} arguments.
}

\usage{
\S4method{overlay}{SpatRaster,SpatRaster}(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{fun}{function to be applied. The first two arguments the function should expect must be numerical vectors}
  \item{...}{additional arguments to be passed to \code{fun}} 
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{ \link[terra]{app},  \link[terra]{math}} }

\examples{
r1  <- rast(ncol=10, nrow=10)
r2  <- rast(ncol=10, nrow=10)
values(r1) <- 1:ncell(r1)
values(r2) <- 10
x <- overlay(r1, r2, fun=function(x,y){return(x/y)})
}
 
\keyword{methods}
\keyword{spatial}
