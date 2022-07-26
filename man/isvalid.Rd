\name{valid}

\alias{is.valid}
\alias{is.valid,SpatVector-method}
\alias{makeValid}
\alias{makeValid,SpatVector-method}

  
\title{Check or fix polygon validity}

\description{
Check the validity of polygons or attempt to fix it
}

\usage{
\S4method{is.valid}{SpatVector}(x, messages=FALSE, as.points=FALSE)

\S4method{makeValid}{SpatVector}(x)
}

\arguments{
\item{x}{SpatVector}
\item{messages}{logical. If \code{TRUE} the error messages are returned}
\item{as.points}{logical. If \code{TRUE}, it is attempted to return locations where polygons are invalid as a SpatVector or points}
}

\value{
logical
}

\examples{
w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
is.valid(w)

w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 4 -2, 0 -5))")
is.valid(w)
is.valid(w, TRUE)

plot(w)
points(cbind(4.54, -2.72), cex=2, col="red")

}
\keyword{spatial}

