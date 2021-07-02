\name{app}

\docType{methods}

\alias{app}
\alias{app,SpatRaster-method}
\alias{app,SpatRasterDataset-method}

\title{Apply a function to the cells of a SpatRaster}

\description{ 
Apply a function to the values of each cell of a SpatRaster. Similar to \code{\link[base]{apply}} -- think of each layer in a SpatRaster as a column (or row) in a matrix. 

This is generally used to summarize the values of multiple layers into one layer; but this is not required.

\code{app} calls function \code{fun} with the raster data as first argument. Depending on the function supplied, the raster data is represented as either a matrix in which each layer is a column, or a vector representing a cell. The function should return a vector or matrix that is divisible by ncell(x). Thus, both "sum" and "rowSums" can be used, but "colSums" cannot be used.

You can also apply a function \code{fun} across datasets by layer of a \code{SpatRasterDataset}. In that case, summarization is across SpatRasters, not across layers.
}

\usage{
\S4method{app}{SpatRaster}(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())

\S4method{app}{SpatRasterDataset}(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{fun}{function that operates on a vector or matrix}
  \item{...}{additional arguments for \code{fun}. These are typically numerical constants. They should *never* be another SpatRaster}  
  \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. Ignored for C++ level implemented functions "max", "min", "mean", "range", "prod", "sum", "any", and "all"}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
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

## SpatRasterDataset
sd <- sds(x, x*2, x/3)
a <- app(sd, max)
a
# same as 
max(x, x*2, x/3)

## also works for a single layer
f <- function(i) (i+1) * 2 * i + sqrt(i)
s <- app(r, f)
# same as above, but that is not memory-safe
# and has no filename argument 
s <- f(r)
}


\keyword{methods}
\keyword{spatial}
