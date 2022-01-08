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


# invalid polygon
\donttest{
w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 4 -2))")
is.valid(w)
}

r <- rast(nrow=10, ncol=10, xmin=0, xmax=1, ymin=0, ymax=1, vals=1)
r[c(44,55)] <- NA
p <- as.polygons(r)
is.valid(p)
# self-intersecting hole
geom(p, TRUE)

x <- is.valid(p, as.points=TRUE)
plot(p, 1)
points(x)
v <- makeValid(p)
is.valid(v)
# the hole has two parts now
geom(v, TRUE)
}

\keyword{spatial}

