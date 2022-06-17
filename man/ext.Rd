\name{ext}

\docType{methods}

\alias{$,SpatExtent-method}
\alias{$<-,SpatExtent-method}
\alias{[,SpatExtent,missing,missing-method}
\alias{[,SpatExtent,numeric,missing-method}
\alias{[<-,SpatExtent,numeric,missing-method}

\alias{ext}
\alias{ext<-}
\alias{ext,SpatExtent-method}
\alias{ext,SpatRaster-method}
\alias{ext,SpatVector-method}
\alias{ext,SpatVectorProxy-method}
\alias{ext,SpatRasterDataset-method}
\alias{ext,SpatRasterCollection-method}
\alias{ext<-,SpatRaster,SpatExtent-method}
\alias{ext<-,SpatRaster,numeric-method}
\alias{ext,numeric-method}
\alias{ext,missing-method}
\alias{ext,Spatial-method}
\alias{ext,Raster-method}
\alias{ext,Extent-method}
\alias{ext,sf-method}


\title{Create, get or set a SpatExtent}

\description{
Get a SpatExtent of a SpatRaster, or coordinates from such an object. Or create a SpatExtent from a vector (length=4; order= xmin, xmax, ymin, ymax)

See \code{\link{set.ext}} to set the extent in place.
}

\usage{
\S4method{ext}{SpatRaster}(x, cells=NULL)

\S4method{ext}{SpatVector}(x)

\S4method{ext}{numeric}(x, ...)

\S4method{ext}{SpatRaster,SpatExtent}(x)<-value

\S4method{ext}{SpatRaster,numeric}(x)<-value

\S4method{$}{SpatExtent}(x, name)

\S4method{$}{SpatExtent}(x, name)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{cells}{postive integer (cell) numbers to subset the extent to area covered by these cells}
  \item{value}{SpatExtent, or numeric vector of lenght four (xmin, xmax, ymin, ymax), or a single number with the \code{$} method}  
  \item{name}{charcter, one of xmin, xmax, ymin, or ymax}  
  \item{...}{if \code{x} is a single numeric value, additional numeric values for xmax, ymin, and ymax}  
}
 
\value{
A \code{\link{SpatExtent}} object. 
}

\examples{
r <- rast()
e <- ext(r)
as.vector(e)
as.character(e)

ext(r) <- c(0, 2.5, 0, 1.5)
r
er <- ext(r)

round(er)
# go "in"
floor(er)
# go "out"
ceiling(er)

ext(r) <- e

}

\keyword{spatial}

