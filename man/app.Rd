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

You can also apply a function \code{fun} across datasets by layer of a \code{SpatRasterDataset}. In that case, summarization is by layer across SpatRasters. 
}

\usage{
\S4method{app}{SpatRaster}(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())

\S4method{app}{SpatRasterDataset}(x, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{fun}{a function that operates on a vector or matrix. This can be a function that is defined in base-R or in a package, or a function you write yourself (see examples). Functions that return complex output (e.g. a list) may need to be wrapped in your own function to simplify the output to a vector or matrix. The following functions have been re-implemented in C++ for speed: "sum", "mean", "median", "modal", "which", "which.min", "which.max", "min", "max", "prod", "any", "all", "sd", "std", "first". To use the base-R function for say, "min", you could use something like \code{fun=function(i) min(i)} or the equivalent \code{fun = \(i) min(i)}}
  \item{...}{additional arguments for \code{fun}. These are typically numerical constants. They should *never* be another SpatRaster}  
  \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object. Ignored for functions that are implemented by terra in C++ (see under fun)}
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

\seealso{ \code{\link{lapp}}, \code{\link{tapp}}, \code{\link[terra]{Math-methods}}, \code{\link{roll}} }


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
x <- c(r, sqrt(r), r+50)
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
# and as (but slower)
b <- app(sd, function(i) max(i))


## also works for a single layer
f <- function(i) (i+1) * 2 * i + sqrt(i)
s <- app(r, f)
# same as above, but that is not memory-safe
# and has no filename argument 
s <- f(r)

\dontrun{
#### multiple cores 
test0 <- app(x, sqrt) 
test1 <- app(x, sqrt, cores=2)

testfun <- function(i) { 2 * sqrt(i) }
test2 <- app(x, fun=testfun, cores =2)

## this fails because testfun is not exported to the nodes
# test3 <- app(x, fun=function(i) testfun(i), cores=2)
## to export it, add it as argument to fun
test3 <- app(x, fun=function(i, ff) ff(i), cores =3, ff=testfun)
}
}


\keyword{methods}
\keyword{spatial}
