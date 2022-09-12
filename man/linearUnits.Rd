\name{linearUnits}

\alias{linearUnits}
\alias{linearUnits,SpatRaster-method}
\alias{linearUnits,SpatVector-method}

  
\title{Linear units of the coordinate reference system}

\description{
Get the linear units of the coordinate reference system (crs) of a SpatRaster or SpatVector expressed in m. The value returned is used internally to transform area and perimeter measures to meters. The value returned for longitude/latitude crs is zero.
}

\usage{
\S4method{linearUnits}{SpatRaster}(x)

\S4method{linearUnits}{SpatVector}(x)

}

\arguments{
  \item{x}{SpatRaster or SpatVector}
}

\value{
numeric (meter)
}

\seealso{\code{\link{crs}}}

\examples{
x <- rast()
crs(x) <- ""
linearUnits(x)

crs(x) <- "+proj=longlat +datum=WGS84"
linearUnits(x)

crs(x) <- "+proj=utm +zone=1 +units=cm"
linearUnits(x)

crs(x) <- "+proj=utm +zone=1 +units=km"
linearUnits(x)

crs(x) <- "+proj=utm +zone=1 +units=us-ft"
linearUnits(x)
}


\keyword{spatial}

