\name{focalReg}

\alias{focalReg}
\alias{focalReg,SpatRaster-method}


\title{Focal regression}

\description{
Calculate the coefficients for a focal ("moving window") OLS regression model. 
}

\usage{
\S4method{focalReg}{SpatRaster}(x, w=3, na.rm=TRUE,  
		fillvalue=NA, weighted=FALSE, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster with at least two layers. The first is the "Y" (dependent) variable and the remainder are the "X" (independent) variables}

\item{w}{numeric or matrix to define the focal window. The window an be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}. Note that if \code{weighted=FALSE} all values that are 0 are set to NA, and all values that are neither 0 or NA are set to 1}

\item{na.rm}{logical. Should missing values be removed?}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{weighted}{logical. If \code{TRUE} and \code{w} is a matrix, the values of w are used as weights. Note that this makes the method much slower}

\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}


\seealso{ \code{\link{focal}}, \code{\link{focalValues}} }


\examples{
r <- rast(ncols=10, nrows=10, ext(0, 10, 0, 10))
values(r) <- 1:ncell(r)
x <- c(r, init(r, runif) * r)
f <- focalReg(x, 3)
}


\keyword{spatial}

