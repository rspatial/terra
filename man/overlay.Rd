\name{overlay}

\docType{methods}

\alias{overlay}
\alias{overlay,SpatRaster,SpatRaster-method}

\title{Overlay SpatRaster objects}

\description{
Create a new SpatRaster object, by combining two or more SpatRaster objects. 

You should supply a function \code{fun} to determine how the SpatRasters are combined. The number of arguments in the function must match the number of SpatRaster objects (or take any number). For example, if you combine two SpatRaster objects you could use multiply: \code{fun=function(x,y){return(x*y)}} percentage: \code{fun=function(x,y){return(100 * x / y)}}. If you combine three layers you could use \code{fun=function(x,y,z){return((x + y) * z)}}

Note that the function must work for vectors (not only for single numbers). That is, it must return the same number of elements as its input vectors. Alternatively, you can also supply a function such as \code{sum}, that takes \code{n} arguments.

}

\usage{
\S4method{overlay}{SpatRaster,SpatRaster}(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{fun}{Function to be applied. The number of arguments of the function should match the number of objects, or it should take any number of arguments. When using multi-layer objects the function should match the number of layers of the object (unless unstack=FALSE) }
  \item{...}{Additional arguments to \code{fun}} 
   \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}

}



\value{
SpatRaster
}

\seealso{\code{ \link[terra]{reduce},  \link[terra]{Arith-methods}} }

\examples{
r1  <- rast(ncol=10, nrow=10)
r2  <- rast(ncol=10, nrow=10)
values(r1) <- 1:ncell(r1)
values(r2) <- ncell(r1):1
x <- overlay(r1, r2, fun=function(x,y){return(x+y)})
}
 
\keyword{methods}
\keyword{spatial}
