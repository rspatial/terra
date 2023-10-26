\name{lapp}

\docType{methods}

\alias{lapp}
\alias{lapp,SpatRaster-method}
\alias{lapp,SpatRasterDataset-method}

\title{Apply a function to layers of a SpatRaster, or sub-datasets of a SpatRasterDataset}

\description{
Apply a function to a SpatRaster, using layers as arguments.

The number of arguments in function \code{fun} must match the number of layers in the SpatRaster (or the number of sub-datasets in the SpatRasterDataset). For example, if you want to multiply two layers, you could use this function: \code{fun=function(x,y){return(x*y)}} percentage: \code{fun=function(x,y){return(100 * x / y)}}. If you combine three layers you could use \code{fun=function(x,y,z){return((x + y) * z)}}

Before you use the function, test it to make sure that it is vectorized. That is, it should work for vectors longer than one, not only for single numbers. Or if the input SpatRaster(s) have multiple layers, it should work for a matrix (multiple cells) of input data (or matrices in the case of a SpatRasterDataSet). The function must return the same number of elements as its input vectors, or multiples of that. Also make sure that the function is \code{NA}-proof: it should returns the same number of values when some or all input values are \code{NA}. And the function must return a vector or a matrix, not a \code{data.frame}. To test it, run it with \code{do.call(fun, data)} (see examples). 

Use \code{\link{app}} for summarize functions such as \code{sum}, that take any number of arguments; and \code{\link{tapp}} to do so for groups of layers. 
}

\usage{
\S4method{lapp}{SpatRaster}(x, fun, ..., usenames=FALSE, cores=1, filename="", overwrite=FALSE, wopt=list())

\S4method{lapp}{SpatRasterDataset}(x, fun, ..., usenames=FALSE, recycle=FALSE, 
    filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{fun}{a function that takes a vector and can be applied to each cell of \code{x}}
  \item{...}{additional arguments to be passed to \code{fun}} 
  \item{usenames}{logical. Use the layer names (or dataset names if \code{x} is a SpatRasterDataset) to match the function arguments? If \code{FALSE}, argument matching is by position}
  \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object}
  \item{recycle}{logical. Recycle layers to match the subdataset with the largest number of layers}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{ \link{app}, \link[terra]{tapp},  \link[terra]{math}} }

\note{
Use \code{\link{sapp}} or \code{lapply} to apply a function that takes a SpatRaster as argument to each layer of a SpatRaster (that is rarely necessary).
}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra")) + 1  
ss <- s[[2:1]]

fvi <- function(x, y){ (x - y ) / (x + y) } 
# test the function
data <- list(c(1:5,NA), 6:1)
do.call(fvi, data)

x <- lapp(ss, fun=fvi )

# which is the same as supplying the layers to "fun"
# in some cases this will be much faster 
y <- fvi(s[[2]], s[[1]])

f2 <- function(x, y, z){ (z - y + 1) / (x + y + 1) } 
p1 <- lapp(s, fun=f2 )

p2 <- lapp(s[[1:2]], f2, z=200)

# the usenames argument

fvi2 <- function(red, green){ (red - green ) / (red + green) } 
names(s)
x1 <- lapp(s[[1:2]], fvi2, usenames=TRUE)
x2 <- lapp(s[[2:1]], fvi2, usenames=TRUE)
# x1 and x2 are the same, despite the change in the order of the layers
# x4 is also the same, but x3 is not
x3 <- lapp(s[[2:1]], fvi2, usenames=FALSE)
x4 <- lapp(s, fvi2, usenames=TRUE)

# while this would fail because 
# there are too many layers in s
# x5 <- lapp(s, fvi2, usenames=FALSE)

pairs(c(x1, x2, x3, x4))

## SpatRasterDataset
x <- sds(s, s[[1]]+50)
fun <- function(x, y) { x/y }

# test "fun"
data <- list(matrix(1:9, ncol=3), matrix(9:1, ncol=3))
do.call(fun, data)

lapp(x, fun, recycle=TRUE)

# the same, more concisely
z <- s / (s[[1]]+50)
}
 
\keyword{methods}
\keyword{spatial}
