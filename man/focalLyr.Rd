\name{focalLyr}

\alias{focalLyr}
\alias{focalLyr,SpatRaster-method}


\title{Focal function across two layers}

\description{
Calculate values for focal regions in two neighboring layers. A function is applied to the first and second layer, then to the second and third layer, etc.
}

\usage{
\S4method{focalLyr}{SpatRaster}(x, w=3, fun, ..., fillvalue=NA, 
		expand=FALSE, filename="", overwrite=FALSE, wopt=list()) 
}

\arguments{
\item{x}{SpatRaster with at least two layers}

\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}}

\item{fun}{A function with at least two arguments (one for each layer)}

\item{...}{additional arguments for \code{fun}}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{expand}{logical. If \code{TRUE} The value of the cells in the virtual rows and columns outside of the raster are set to be the same as the value on the border}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\seealso{ \code{\link{layerCor}}, \code{\link{focal}}, \code{\link{focalValues}} }


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))
set.seed(0)
r[[1]] <- flip(r[[1]], "horizontal")
r[[2]] <- flip(r[[2]], "vertical") + init(rast(r,1), runif)
r[[3]] <- init(rast(r,1), runif)
x <- focalLyr(r, w=5, cor)

x <- focalLyr(r, w=3, function(x, y) suppressWarnings(cor(x, y)))
plot(x)
}



\keyword{spatial}

