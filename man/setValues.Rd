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

\S4method{setValues}{SpatRaster,ANY}(x, values, keeptime=TRUE, keepunits=TRUE, keepnames=FALSE, props=FALSE)

\S4method{values}{SpatVector,ANY}(x)<-value
}


\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{value}{For SpatRaster: numeric, matrix or data.frame. The length of the numeric values must match the total number of cells (ncell(x) * nlyr(x)), or be a single value. The number of columns of the matrix or data.frame must match the number of layers of \code{x}, and the number of rows must match the number of cells of \code{x}. It is also possible to use a matrix with the same number of rows as \code{x} and the number of columns that matches \code{ncol(x) * nlyr(x)}. 

 For SpatVector: data.frame, matrix, vector, or NULL}
 \item{values}{Same as for \code{value}}
 \item{keeptime}{logical. If \code{TRUE} the time stamps are kept}
 \item{keepunits}{logical. If \code{FALSE} the units are discarded}
 \item{keepnames}{logical. If \code{FALSE} the layer names are replaced by the column names in \code{y} (if present)}
 \item{props}{logical. If \code{TRUE} the properties (categories and color-table) are kept}
}


\value{
The same object type as \code{x}
}


\seealso{\code{\link{values}}, \code{\link{init}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- setValues(r, 1:ncell(r))
x
values(x) <- runif(ncell(x))
x
head(x)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
values(v) <- data.frame(ID=1:12, name=letters[1:12])
head(v)
}

\keyword{spatial}
\keyword{methods}
