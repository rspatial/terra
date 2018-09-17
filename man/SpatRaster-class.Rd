\name{SpatRaster-class}

\docType{class}

\alias{SpatRaster}
\alias{SpatRaster-class}
\alias{Rcpp_SpatRaster-class}

\alias{RasterSource}
\alias{RasterSource-class}
\alias{Rcpp_RasterSource-class}


\alias{show,SpatRaster-method}


\title{ SpatRaster class}

\description{
A raster is a database organized as a rectangular grid that is sub-divided into rectangular cells of equal area (in terms of the units of the coordinate reference system) and layers. 

An object of the \code{SpatRaster} class can point to one or more files on disk that holds the values of the raster cells, or hold these values in memory. Or it can not have any associated values at all.

Objects can be created with the \code{\link{rast}} function.
 
This package defines the SpatRaster "S4 class" to manipulate such data. The R object only contains a pointer to the C++ class "Rcpp_SpatRaster". 
}


\section{Objects from the Class}{
Objects can be created with the \code{rast} function.
}

\section{Slots}{
Slots for SpatRaster objects
	\describe{
    \item{\code{ptr}:}{pointer to the C++ class} 
}
}


\examples{
rast()
}

\keyword{classes}
\keyword{spatial}
