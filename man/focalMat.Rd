\name{focalMat}

\alias{focalMat}

\title{Focal weights matrix}

\description{
Make a focal ("moving window") weight matrix for use in the \code{\link{focal}} function. The sum of the values adds up to one.
}

\usage{
focalMat(x, d, type=c('circle', 'Gauss', 'rectangle'), fillNA=FALSE)
}

\arguments{
\item{x}{SpatRaster}
\item{d}{numeric. If \code{type=circle}, the radius of the circle (in units of the crs). If \code{type=rectangle} the dimension of the rectangle (one or two numbers). If \code{type=Gauss} the size of sigma, and optionally another number to determine the size of the matrix returned (default is 3*sigma)} 
\item{type}{character indicating the type of filter to be returned}
\item{fillNA}{logical. If \code{TRUE}, zeros are set to \code{NA} such that they are ignored in the computations. Only applies to \code{type="circle"}}
}

\value{
matrix that can be used with \code{\link{focal}}
}


\examples{
r <- rast(ncols=180, nrows=180, xmin=0)
focalMat(r, 2, "circle")

focalMat(r, c(2,3), "rect")

# Gaussian filter for square cells
gf <- focalMat(r, 1, "Gauss")
}


\keyword{spatial}

