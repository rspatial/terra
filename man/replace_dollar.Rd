\name{replace_dollar}

\alias{$<-}
\alias{$<-,SpatExtent-method}
\alias{$<-,SpatVector-method}
\alias{$<-,SpatRaster-method}


\title{Replace with $<-}

\description{
Replace a layer of a SpatRaster, or an attribute variable of a SpatVector
}

\usage{
\S4method{$}{SpatRaster}(x, name) <- value

\S4method{$}{SpatVector}(x, name)<-value

\S4method{$}{SpatExtent}(x, name) <- value
}


\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{name}{character. If \code{x} is a SpatRaster: layer name. If \code{x} is a SpatVector: variable name. If \code{x} is a SpatExtent: "xmin", "xmax". "ymin" or "ymax"}
  \item{value}{if \code{x} is a SpatRaster, a SpatRaster for which this \code{TRUE}: \code{nlyr(value) == length(i)}; if \code{x} is a SpatVector, a vector of new values; if \code{x} is a SpatExtent a single number}
}


\value{
Same as \code{x}
}


\seealso{
\code{ \link{[[<-}, \link{[<-}, \link{$}}
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$ID_1 <- LETTERS[1:12]
v$new <- sample(12)
values(v)
}

\keyword{ spatial }
