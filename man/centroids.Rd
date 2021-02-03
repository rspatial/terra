\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}


\title{Get centroids}

\description{
Get the centroids for the polygons of a SpatVector
}

\usage{
\S4method{centroids}{SpatVector}(x)
}

\arguments{
\item{x}{SpatVector}
}


\value{SpatVector of points}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
}

\keyword{spatial}
