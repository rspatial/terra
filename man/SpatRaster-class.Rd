\name{SpatRaster-class}

\docType{class}

\alias{SpatRaster}
\alias{SpatRaster-class}
\alias{Rcpp_SpatRaster-class}
\alias{PackedSpatRaster-class}

\alias{RasterSource}
\alias{RasterSource-class}
\alias{Rcpp_RasterSource-class}

\alias{SpatCategories}
\alias{SpatCategories-class}
\alias{Rcpp_SpatCategories-class}

\alias{show,SpatRaster-method}


\title{ SpatRaster class}

\description{
A SpatRaster represents a rectangular part of the world that is sub-divided into rectangular cells of equal area (in terms of the units of the coordinate reference system). For each cell can have multiple values ("layers"). 

An object of the \code{SpatRaster} class can point to one or more files on disk that hold the cell values, and/or it can hold these values in memory. These objects can be created with the \code{\link{rast}} method. 

The underlying C++ class "Rcpp_SpatRaster" is not intended for end-users. It is for internal use within this package only.
}


\examples{
rast()
}

\keyword{classes}
\keyword{spatial}
