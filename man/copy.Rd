\name{copy}

\docType{methods}

\alias{copy}
\alias{copy,SpatRaster-method}
\alias{copy,SpatVector-method}


\title{Deep copy}

\description{ 
Make a deep copy of a SpatRaster. This is occaisionaly useful when wanting to use a replacement function in a shallow copy. If use a replacement function to change an object, its shallow copies will also be changed.
}

\usage{
\S4method{copy}{SpatRaster(x, ...)}
}


\arguments{
  \item{x}{SpatRaster}
 \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
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

y <- copy(r)
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
