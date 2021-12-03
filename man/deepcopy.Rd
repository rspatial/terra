\name{deepcopy}

\docType{methods}

\alias{deepcopy}
\alias{deepcopy,SpatRaster-method}
\alias{deepcopy,SpatVector-method}


\title{Deep copy}

\description{ 
Make a deep copy of a SpatRaster or SpatVector. This is occasionally useful when wanting to use a replacement function in a shallow copy. That is a copy that was created like this: \code{x <- y}. If you use a replacement function to change an object, its shallow copies also change.
}

\usage{
\S4method{deepcopy}{SpatRaster}(x)

\S4method{deepcopy}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
}

\value{
Same as \code{x}
}


\examples{
r <- rast(ncols=10, nrows=10, nl=3)
tm <- as.Date("2001-05-03") + 1:3
time(r) <- tm
time(r)
x <- r
time(x) <- tm + 365
time(x)
time(r)

y <- deepcopy(r)
time(y) <- tm - 365
time(y)
time(r)

# or make a new object like this
z <- rast(r)
time(z) <- tm
time(z)
time(r)
}


\keyword{methods}
\keyword{spatial}
