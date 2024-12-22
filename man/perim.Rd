\name{perim}

\alias{perim}
\alias{perim,SpatVector-method}

\alias{perimeter}
\alias{perimeter,SpatVector-method}
  
\title{Perimeter or length}

\description{
This method returns the length of lines or the perimeter of polygons.

When the coordinate reference system is not longitude/latitude, you may get more accurate results by first transforming the data to longitude/latitude with \code{\link{project}}
}

\usage{
\S4method{perim}{SpatVector}(x)
}


\arguments{
  \item{x}{SpatVector}
}

\value{
numeric (m)
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
perim(v)
}

\keyword{methods}
\keyword{spatial}
