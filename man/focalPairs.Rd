\name{focalPairs}

\alias{focalCor}
\alias{focalCor,SpatRaster-method}
\alias{focalPairs}
\alias{focalPairs,SpatRaster-method}


\title{Focal function across two layers}

\description{
Calculate values such as a correlation coefficient for focal regions in two neighboring layers. A function is applied to the first and second layer, then to the second and third layer, etc.
}

\usage{
\S4method{focalPairs}{SpatRaster}(x, w=3, fun, ..., fillvalue=NA, 
		filename="", overwrite=FALSE, wopt=list()) 
}

\arguments{
\item{x}{SpatRaster with at least two layers}

\item{w}{numeric or matrix to define the focal window. The window an be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}. Note that if a matrix with numbers other than zero or one are used, the values are used as weights. For this to work, \code{fun} must have an argument \code{weights}}

\item{fun}{a function with at least two arguments (one for each layer). There is a built-in function "pearson" (for both the weighted and the unweighted Pearson correlation coefficient. This function has an additional argument \code{na.rm=FALSE}}

\item{...}{additional arguments for \code{fun}}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\seealso{ \code{\link{layerCor}}, \code{\link{focalReg}}, \code{\link{focal}}, \code{\link{focal3D}}}


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))
set.seed(0)
r[[1]] <- flip(r[[1]], "horizontal")
r[[2]] <- flip(r[[2]], "vertical") + init(rast(r,1), runif)
r[[3]] <- init(rast(r,1), runif)

x <- focalPairs(r, w=5, "pearson", na.rm=TRUE)
plot(x)


# suppress warning "the standard deviation is zero"
suppressWarnings(x <- focalPairs(r, w=5, "pearson", use="complete.obs"))

z <- focalPairs(r, w=9, function(x, y) mean(x) + mean(y))
}



\keyword{spatial}

