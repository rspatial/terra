\name{na.omit}

\alias{na.omit}
\alias{na.omit,SpatVector-method} 

\alias{is.na,SpatVector-method} 
 
 
\title{Find and remove geometries that are NA}

\description{
Find geometries that are NA; or remove geometries and/or records that are \code{NA}. 
}


\usage{
\S4method{is.na}{SpatVector}(x)

\S4method{na.omit}{SpatVector}(object, field=NA, geom=FALSE)
}

\arguments{
  \item{x}{SpatVector}  
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
