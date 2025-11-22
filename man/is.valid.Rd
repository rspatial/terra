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

\S4method{makeValid}{SpatVector}(x, buffer=FALSE)

\S4method{is.valid}{SpatExtent}(x)
}

\arguments{
\item{x}{SpatVector or SpatExtent}
\item{messages}{logical. If \code{TRUE} the error messages are returned}
\item{as.points}{logical. If \code{TRUE}, it is attempted to return locations where polygons are invalid as a SpatVector or points}
\item{buffer}{logical. If \code{TRUE} the zero-width buffer method is used to create valid polygons. Be careful when  using this method because it may result in data loss. For example, only a single part of a self-intersecting may be preserved. See the example below}
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
v1 <- makeValid(w)
is.valid(v1)
v2 <- makeValid(w, buffer=TRUE)

plot(disagg(v1), col=c("light blue", "gray"))
# valid but incomplete
lines(v2, col="red", lwd=3)

plot(w)
points(cbind(4.54, -2.72), cex=2, col="red")

e <- ext(0, 1, 0, 1)
is.valid(e)

ee <- ext(0, 0, 0, 0)
is.valid(ee)

}
\keyword{spatial}

