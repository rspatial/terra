\name{is.valid}

\alias{is.valid}
\alias{is.valid,SpatVector-method}

  
\title{Check polygon validity}

\description{
Check the validity of polygons
}

\usage{
\S4method{is.valid}{SpatVector}(x, messages=FALSE, as.points=FALSE)
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

# invalid polygon
#w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 4 -2))")
#is.valid(w)

}

\keyword{spatial}

