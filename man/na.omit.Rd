\name{na.omit}

\alias{na.omit}
\alias{na.omit,SpatVector-method} 
 
 
\title{na.omit for SpatVector}

\description{
Remove empty geometries and/or records that are \code{NA} from a SpatVector. 
}


\usage{
\S4method{na.omit}{SpatVector}(object, field=NA, geom=FALSE)
}

\arguments{
  \item{object}{SpatVector}  
  \item{field}{character or NA. If \code{NA}, missing values in the attributes are ignored. Other values are either one or more field (variable) names, or \code{""} to consider all fields}
  \item{geom}{logical. If \code{TRUE} empty geometries are removed}  
}
 
\value{
SpatVector
}
 
\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$test <- c(1,2,NA)
nrow(v)
x <- na.omit(v, "test")
nrow(x)
}

\keyword{spatial}
