\name{focalReg}

\alias{focalReg}
\alias{focalReg,SpatRaster-method}

\title{Focal regression}

\description{
Calculate values for a moving-window by comparing the value in one layers with the values in one to many other layers. A typical case is the computation of the coefficients for a focal linear regression model. 
}

\usage{
\S4method{focalReg}{SpatRaster}(x, w=3, fun="ols", ..., fillvalue=NA, filename="", overwrite=FALSE, wopt=list()) 
}

\arguments{
\item{x}{SpatRaster with at least two layers. The first is the "Y" (dependent) variable and the remainder are the "X" (independent) variables}

\item{w}{numeric or matrix to define the focal window. The window an be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}. Note that if a matrix with numbers other than zero or one are used, the values are used as weights. For this to work, \code{fun} must have an argument \code{weights}}

\item{fun}{a function with at least two arguments (one for each layer). There is a built-in function "ols" for both the weighted and unweighted Ordinary Least Square regression. This function has an additional argument \code{na.rm=FALSE} and \code{intercept=TRUE}}

\item{...}{additional arguments for \code{fun}}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{wopt}{additional arguments for writing files as in \code{\link{writeRaster}}}

}

\value{SpatRaster}


\seealso{ \code{\link{focal}, \link{focal3D}}, \link{focalValues}} 


\examples{
r <- rast(ncols=10, nrows=10, ext(0, 10, 0, 10))
values(r) <- 1:ncell(r)
x <- c(r, init(r, runif) * r)
f <- focalReg(x, 3)
}


\keyword{spatial}

