\name{datatype}

\docType{methods}

\alias{datatype}
\alias{datatype,SpatVector-method}
\alias{datatype,SpatRaster-method}

\title{Data type of a SpatRaster or SpatVector}

\description{
Get the data types of the fields (attributes, variables) of a SpatVector or of the file(s) associated with a SpatRaster. A (layer of a) SpatRaster has no datatype if it has no values, or if the values are in memory. 
} 

\usage{
\S4method{datatype}{SpatRaster}(x, bylyr=TRUE)

\S4method{datatype}{SpatVector}(x)
}

\arguments{
\item{x}{SpatRaster or SpatVector}

\item{bylyr}{logical. If \code{TRUE} a value is returned for each layer. Otherwise, a value is returned for each data source (such as a file)}
}

\value{character}

\seealso{
 \code{\link[=is.bool]{Raster data types}} to check / set the type of SpatRaster values.
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
datatype(v)

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
datatype(r)

# no data type
datatype(rast()) 
}

\keyword{methods}
\keyword{spatial}

