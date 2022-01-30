\name{focalCor}

\alias{focalCor}
\alias{focalCor,SpatRaster-method}


\title{Focal function across two layers}

\description{
Calculate values such as a correlation coefficient for focal regions in two neighboring layers. A function is applied to the first and second layer, then to the second and third layer, etc.
}

\usage{
\S4method{focalCor}{SpatRaster}(x, w=3, fun, ..., fillvalue=NA, 
		filename="", overwrite=FALSE, wopt=list()) 
}

\arguments{
\item{x}{SpatRaster with at least two layers}

\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}}

\item{fun}{A function with at least two arguments (one for each layer)}

\item{...}{additional arguments for \code{fun}}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\seealso{ \code{\link{layerCor}}, \code{\link{focalReg}}, \code{\link{focal}}}


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))
set.seed(0)
r[[1]] <- flip(r[[1]], "horizontal")
r[[2]] <- flip(r[[2]], "vertical") + init(rast(r,1), runif)
r[[3]] <- init(rast(r,1), runif)

# suppress warning "the standard deviation is zero"
suppressWarnings(x <- focalCor(r, w=5, cor))

# this warning does not occur when using a larger window
x <- focalCor(r, w=9, function(x, y) cor(x, y))
plot(x)
}



\keyword{spatial}

