\name{replace_layers}

\docType{methods}

\alias{[[<-}

\alias{[[<-,SpatRaster,character-method}
\alias{[[<-,SpatRaster,numeric-method}
\alias{[[<-,SpatVector,character-method}
\alias{[[<-,SpatVector,numeric-method}


\title{Replace layers or variables}

\description{
Replace the layers of SpatRaster with (layers from) another SpatRaster or replace variables of a SpatVector. You can also create new layers/variables with these methods.
}

\value{SpatRaster}

\usage{
\S4method{[[}{SpatRaster,numeric}(x, i) <- value

\S4method{[[}{SpatRaster,character}(x, i) <- value

\S4method{[[}{SpatVector,numeric}(x, i) <- value

\S4method{[[}{SpatVector,character}(x, i) <- value
}

\arguments{
  \item{x}{SpatRaster or SpatVector}

  \item{i}{if \code{x} is a SpatRaster: layer number(s) of name(s). If \code{x} is a SpatVector: variable  number(s) or name(s) (column of the attributes)}
  
  \item{value}{if \code{x} is a SpatRaster: SpatRaster for which this \code{TRUE}: \code{nlyr(value) == length(i)}. if \code{x} is a SpatVector: vector or data.frame}
}


\seealso{\code{\link{$<-}, \link{[<-}}}

\examples{
# raster
s <- rast(system.file("ex/logo.tif", package="terra"))   
s[["red"]] <- mean(s)
s[[2]] <- sqrt(s[[1]])

# vector
v <- vect(system.file("ex/lux.shp", package="terra"))
v[["ID_1"]] <- 12:1
}

\keyword{methods}
\keyword{spatial}

