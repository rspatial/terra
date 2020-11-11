\name{app}

\docType{methods}

\alias{app}
\alias{app,SpatRaster-method}
\alias{app,SpatDataSet-method}

\title{Apply a function to the cells of a SpatRaster}

\description{ 
Apply a function to values of each cell of a SpatRaster. Similar to \code{\link[base]{apply}} -- think of each layer in a SpatRaster as a column (or row) in a matrix. This is generally used to summarize the values of multiple layers into one layer; but this is not required.

You can also apply a function \code{fun} across datasets by layer of a \code{SpatDataSet}. In that case, summarization is across SpatRasters, not across layers.
}

\usage{
\S4method{app}{SpatRaster}(x, fun, ..., nodes=1, filename="", overwrite=FALSE, wopt=list())

\S4method{app}{SpatDataSet}(x, fun, ..., nodes=1, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatDataSet}
  \item{fun}{function}
  \item{...}{additional arguments for \code{fun}}  
  \item{nodes}{positive integer. If \code{nodes > 1}, a 'parallel' package cluster with that many nodes is created. Ignored for C++ level implemented functions "max", "min", "mean", "range", "prod", "sum", "any", and "all"}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\details{
To speed things up, parallelization is supported, but this is often not helpful, and it may actually be slower. There is only a speed gain if you have many cores (> 8) and/or a very complex (slow) function \code{fun}. If you write \code{fun} yourself, consider supplying a \code{cppFunction} made with the Rcpp package instead (or go have a cup of tea while the computer works for you).
}

\seealso{ \code{\link{lapp}}, \code{\link{tapp}}, \code{\link[terra]{math}} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r, sqrt(r), r-50)
s <- app(x, fun=sum)
s
# for a few generic functions like 
# "sum", "mean", and "max" you can also do
sum(x)

## SpatDataSet
sd <- sds(x, x*2, x/3)
a <- app(sd, max)
a
# same as 
max(x, x*2, x/3)

}


\keyword{methods}
\keyword{spatial}
