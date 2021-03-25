\name{setValues}

\docType{methods}

\alias{values<-}
\alias{values<-,SpatRaster,ANY-method}

\alias{setValues}
\alias{setValues,SpatRaster-method}
\alias{setValues,SpatRaster,ANY-method}

\alias{values<-,SpatVector,data.frame-method}
\alias{values<-,SpatVector,matrix-method}
\alias{values<-,SpatVector,ANY-method}
\alias{values<-,SpatVector,NULL-method}

\alias{setValues,SpatVector-method}
\alias{setValues,SpatVector,ANY-method}

  
\title{Set the values of raster cells or of geometry attributes}

\description{
Set cell values of a SpatRaster or the attributes of a SpatVector. For large SpatRaster objects use \code{\link{init}} instead to set values. 
}

\usage{
\S4method{values}{SpatRaster,ANY}(x)<-value

\S4method{setValues}{SpatRaster,ANY}(x, values)

\S4method{values}{SpatVector,ANY}(x)<-value
}


\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{value}{For SpatRaster: matrix or numeric, the length must match the total number of cells (ncell(x) * nlyr(x)), or be a single value. For SpatVector: data.frame, matrix, vector, or NULL}
 \item{values}{Same as for \code{value}}
}


\value{
setValues returns the same object type as \code{x}
}


\seealso{\code{\link{values}}, \code{\link{init}}}

\examples{
r <- rast(system.file("ex/test.tif", package="terra"))
x <- setValues(r, 1:ncell(r))
values(x) <- runif(ncell(x))
head(x)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
values(v) <- data.frame(ID=1:12, name=letters[1:12])
head(v)
}

\keyword{spatial}
\keyword{methods}
