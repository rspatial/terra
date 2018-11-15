\name{reduce}

\docType{methods}

\alias{reduce}
\alias{reduce,SpatRaster,function-method}

\title{reduce}

\description{ 
reduce values of a SpatRaster, using a formula. 
}

\usage{
\S4method{reduce}{SpatRaster,function}(x, fun, ..., filename="", overwrite=FALSE, writeopt)
}

\arguments{
  \item{x}{SpatRaster object}
  \item{fun}{function}
  \item{...}{Additional arguments for \code{fun}}  
   \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{writeopt}{list. Options for writing files as in \code{\link{writeRaster}}}

}

\value{
SpatRaster
}

\seealso{ \code{ \link[terra]{overlay}}, \link[terra]{Arith-methods}, \link[terra]{Math-methods}}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r,r,r)
s <- reduce(x, fun=sum)
}


\keyword{methods}
\keyword{spatial}
