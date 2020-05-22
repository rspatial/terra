\name{app}

\docType{methods}

\alias{app}
\alias{app,SpatRaster-method}

\title{Apply a function to the cells of a SpatRaster}

\description{ 
Apply a function to values of each cell of a SpatRaster. This is similar to \code{\link[base]{apply}} -- think of each layer in a SpatRaster as a column (or row) in a matrix.

To speed things up, you may use parallel computing, but this is often not helpful, and it may actually be slower. There is only a speed gain if you have many cores (> 8) and/or a very complex (slow) function \code{fun}. If you write \code{fun} yourself, consider supplying a \code{cppFunction} made with the Rcpp package instead (or go have a cup of tea while the computer works for you).
}

\usage{
\S4method{app}{SpatRaster}(x, fun, ..., cls=NULL, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster object}
  \item{fun}{function}
  \item{...}{additional arguments for \code{fun}}  
  \item{cls}{a cluster object created by the 'parallel' package. Ignored for C++ level implemented functions "max", "min", "mean", "range", "prod", "sum", "any", and "all"}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{ \link[terra]{overlay}}, \link[terra]{math} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r, r, r)
s <- app(x, fun=sum)

library(parallel)
detectCores()
cl <- makeCluster(2)
ss <- app(x, fun=sum, cls=cl)


}


\keyword{methods}
\keyword{spatial}
