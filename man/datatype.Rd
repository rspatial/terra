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


\details{
Setting the data type is useful if you want to write values to disk with \code{\link{writeRaster}}. In other cases you can use functions such as \code{round} and \code{floor}, or \code{\link{as.bool}}

raster datatypes are described by 5 characters. The first three indicate whether the values are integer or decimal values. The fourth character indicates the number of bytes used to save the values on disk, and the last character indicates whether the numbers are signed (that is, can be negative and positive values) or not (only zero and positive values allowed)

The following raster datatypes are available:

\tabular{lll}{
\bold{Datatype definition} \tab \bold{minimum possible value} \tab \bold{maximum possible value} \cr
\code{INT1U} \tab 0 \tab  255 \cr
\code{INT2U} \tab 0 \tab  65,534 \cr
\code{INT4U} \tab 0 \tab 4,294,967,296 \cr
\code{INT8U} \tab 0 \tab  18,446,744,073,709,551,616\cr
\code{INT2S} \tab -32,767\tab  32,767 \cr
\code{INT4S} \tab -2,147,483,647 \tab 2,147,483,647 \cr
\code{INT8S} \tab -9,223,372,036,854,775,808 \tab 9,223,372,036,854,775,808 \cr
\code{FLT4S} \tab -3.4e+38 \tab  3.4e+38 \cr
\code{FLT8S} \tab -1.7e+308 \tab   1.7e+308 \cr
}

For all integer types, except the single byte types, the lowest (signed) or highest (unsigned) value is used to store \code{NA}. 

Note that very large integer numbers may be inpreceise as they are internally represented as decimal numbers.

\code{INT4U} is available but they are best avoided as R does not support 32-bit unsigned integers.
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

