\name{is.empty}

\alias{is.empty}
\alias{is.empty,SpatVector-method}
\alias{is.empty,SpatExtent-method}

  
\title{Check if a SpatExtent or SpatVector is empty}

\description{
An empty SpatExtent has no area

An empty SpatVector has no geometries.
}

\usage{
\S4method{is.empty}{SpatExtent}(x)

\S4method{is.empty}{SpatVector}(x)
}

\arguments{
\item{x}{SpatVector or SpatExtent}
}

\value{
logical
}

\examples{
e <- ext(0,0,0,0)
is.valid(e)
is.empty(e)

v <- vect()
is.valid(v)
is.empty(v)
}
\keyword{spatial}

