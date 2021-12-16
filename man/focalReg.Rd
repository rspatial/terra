\name{focalReg}

\alias{focalReg}
\alias{focalReg,SpatRaster-method}


\title{Focal regression}

\description{
Calculate the coefficients for a focal ("moving window") OLS regression model. 
}

\usage{
\S4method{focalReg}{SpatRaster}(x, w=3, na.rm=TRUE,  
		fillvalue=NA, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster with at least two layers. The first is the "Y" (dependent) variable and the remainder are the "X" (independent) variables}

\item{w}{window. The window can be defined as one (for a square) or two numbers (row, col); or with an odd-sized weights matrix. See the Details section in \code{\link{focal}}}

\item{na.rm}{logical. Should missing values be removed?}

\item{fillvalue}{numeric. The value of the cells in the virtual rows and columns outside of the raster}

\item{expand}{logical. If \code{TRUE} The value of the cells in the virtual rows and columns outside of the raster are set to be the same as the value on the border}

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

