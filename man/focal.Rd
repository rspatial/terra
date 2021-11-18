\name{focal}

\alias{focal}
\alias{focal,SpatRaster-method}


\title{Focal values}

\description{
Calculate focal ("moving window") values for the neighborhood of focal cells using a matrix of weights, perhaps in combination with a function.
}

\usage{
\S4method{focal}{SpatRaster}(x, w=3, fun="sum", na.rm=TRUE, na.only=FALSE, 
		fillvalue=NA, expand=FALSE, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster}
  
\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See Details.}

\item{fun}{function that takes multiple numbers, and returns a numeric vector (one or multiple numbers). For example mean, modal, min or max. It should also accept a \code{na.rm} argument, either as actual argument or through use of \code{...}}

\item{na.rm}{logical. Should missing values be removed?}

\item{na.only}{logical. Should only missing values in \code{x} be changed?}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{expand}{logical. If \code{TRUE} The value of the cells in the virtual rows and columns outside of the raster are set to be the same as the value on the border}

\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\details{
\code{focal} 

The window used must have odd dimensions. If you need even sides, you can use a matrix and add a column or row with weights of zero.

Window values are typically 0 or 1, or a value between 0 and 1 if you are using a rectangular area and/or the "sum" function. They can also be \code{NA}; these are ignored in the computation. That can be useful to compute the mean, min, or max value for a non-rectangular area.  

Example weight matrices

Laplacian filter: \code{filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)}

Sobel filters (for edge detection):
\code{fx=matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3)}
\code{fy=matrix(c(1,0,-1,2,0,-2,1,0,-1), nrow=3)}
}


\value{SpatRaster}


\seealso{ \code{\link{focalMat}}, \code{\link{focalValues}} }


\examples{
r <- rast(ncols=10, nrows=10, ext(0, 10, 0, 10))
values(r) <- 1:ncell(r)

f <- focal(r, w=3, fun=function(x, ...) quantile(x, c(.25, .5, .75), ...), na.rm=TRUE) 

f <- focal(r, w=3, fun="mean") 

# the following two statements are equivalent:
a <- focal(r, w=matrix(1/9, nc=3, nr=3))
b <- focal(r, w=3, fun=mean, na.rm=FALSE)

# but this is different
d <- focal(r, w=3, fun=mean, na.rm=TRUE)
}


\keyword{spatial}

