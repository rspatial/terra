\name{replace_layers}

\docType{methods}

\alias{$<-,SpatRaster-method}
\alias{[[<-,SpatRaster,character,missing-method}
\alias{[[<-,SpatRaster,numeric,missing-method}

\title{Replace layers of a SpatRaster}

\description{
Replace the layers of SpatRaster with (layers from) another SpatRaster
}

\value{SpatRaster}

\usage{
\S4method{[[}{SpatRaster,numeric}(x, i) <- value

\S4method{[[}{SpatRaster,character}(x, i) <- value

\S4method{$}{SpatRaster,character}(x, name) <- value
}

\arguments{
  \item{x}{SpatRaster}
  \item{i}{layer number(s) of name(s)}
  \item{name}{layer name}
  \item{value}{SpatRaster for which this \code{TRUE}: \code{nlyr(value) == length(i)}}
}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
s[["red"]] <- mean(s)
s[[2]] <- sqrt(s[[1]])
}

\keyword{methods}
\keyword{spatial}

