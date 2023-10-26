\name{is.valid}

\alias{is.valid}
\alias{is.valid,SpatVector-method}
\alias{is.valid,SpatExtent-method}
\alias{makeValid}
\alias{makeValid,SpatVector-method}

  
\title{Check or fix polygon or extent validity}

\description{
Check the validity of polygons or attempt to fix it. Or check the validity of a SpatExtent.
}

\usage{
\S4method{is.valid}{SpatVector}(x, messages=FALSE, as.points=FALSE)

\S4method{makeValid}{SpatVector}(x)

\S4method{is.valid}{SpatExtent}(x)
}

\arguments{
\item{x}{SpatVector or SpatExtent}
\item{messages}{logical. If \code{TRUE} the error messages are returned}
\item{as.points}{logical. If \code{TRUE}, it is attempted to return locations where polygons are invalid as a SpatVector or points}
}

\value{
logical
}

\seealso{\code{\link{topology}}}

\examples{
w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
is.valid(w)

w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 4 -2, 0 -5))")
is.valid(w)
is.valid(w, TRUE)

plot(w)
points(cbind(4.54, -2.72), cex=2, col="red")


e <- ext(0, 1, 0, 1)
is.valid(e)

ee <- ext(0, 0, 0, 0)
is.valid(ee)

}
\keyword{spatial}

