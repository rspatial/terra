\name{deepcopy}

\docType{methods}

\alias{deepcopy}
\alias{deepcopy,SpatRaster-method}
\alias{deepcopy,SpatVector-method}


\title{Deep copy}

\description{ 
Make a deep copy of a SpatRaster or SpatVector. This is occasionally useful when using an in=place replacement function that does not make copy, such as \code{\link{set.ext}}. 
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


?\examples{
r <- rast(ncols=10, nrows=10, nl=3)
x <- r
y <- deepcopy(r)
ext(r)
set.ext(x, c(0,10,0,10))
ext(x)
ext(r)
ext(y)
}


\keyword{methods}
\keyword{spatial}
