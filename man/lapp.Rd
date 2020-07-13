\name{lapp}

\docType{methods}

\alias{lapp}
\alias{lapp,SpatRaster-method}
\alias{lapp,SpatDataSet-method}

\title{Apply a function to layers of a SpatRaster, or to sub-datasets of a SpatDataSet}

\description{
Apply a function to layers of a SpatRaster ("overlay"). 

The number of arguments in function \code{fun} must match the number of layers in the SpatRaster (or the number of sub-datasets in the SpatDataSet). For example, if you want to multiply two layers, you could use this function : \code{fun=function(x,y){return(x*y)}} percentage: \code{fun=function(x,y){return(100 * x / y)}}. If you combine three layers you could use \code{fun=function(x,y,z){return((x + y) * z)}}

Before you use the function, test it to make sure that it works for vectors (not only for single numbers). That is, it must return the same number of elements as its input vectors, or multiples of that. Make sure it also returns the same number of values when some or all input values are \code{NA}. Also, the function must return a vector or matrix, not a \code{data.frame}.

Use \code{\link{app}} for summarize functions such as \code{sum}, that take any number of arguments; and \code{\link{tapp}} to so for groups of layers. 

}

\usage{
\S4method{lapp}{SpatRaster}(x, fun, ..., usenames=FALSE, filename="", overwrite=FALSE, wopt=list())

\S4method{lapp}{SpatDataSet}(x, fun, ..., recycle=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatDataSet}
  \item{fun}{function to be applied}
  \item{...}{additional arguments to be passed to \code{fun}} 
  \item{usenames}{logical. Use the layer names to match the function arguments? If \code{FALSE} matching is by position}
  \item{recycle}{logical. Recycle layers to match the subdataset with the largest number of layers}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{ \link[terra]{app}, \link[terra]{tapp},  \link[terra]{math}} }

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
ss <- s[[2:1]]

fvi <- function(x, y){ (x - y ) / (x + y) } 
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

# while this would give an error because 
# there are too many layers in s
# x5 <- lapp(s, fvi2, usenames=FALSE)

pairs(c(x1, x2, x3, x4))

## SpatDataSet
x <- sds(s, s[[1]]+50)
lapp(x, function(x, y) x/y, recycle=TRUE)
}
 
\keyword{methods}
\keyword{spatial}
