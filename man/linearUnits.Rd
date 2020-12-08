\name{linearUnits}

\alias{linearUnits}
\alias{linearUnits,SpatRaster-method}
\alias{linearUnits,SpatVector-method}

  
\title{Linear units of the coordinate reference system}

\description{
Get the linear units of the coordinate reference system (crs) of a SpatRaster or SpatVector. 
The value returned can be used to multiply linear distances to transform them to meters.
The value returned for longitude/latitude crs is zero.
}

\usage{
\S4method{linearUnits}{SpatRaster}(x, ...)

\S4method{linearUnits}{SpatVector}(x, ...)

}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{...}{additional arguments. None implemented}
}

\value{
numeric  
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
}


\keyword{spatial}

