\name{perim}

\alias{perim}
\alias{perim,SpatVector-method}

\alias{perimeter}
\alias{perimeter,SpatVector-method}
  
\title{Perimeter or length}

\description{
This method returns the length of lines or the perimeter of polygons.

When the crs is not longitude/latitude, you may get more accurate results by first unprojecting the SpatVector (you can use \code{\link{project}} to transform the crs to longitude/latitude)
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
