\name{GeoExtent-class}

\docType{class}

\alias{GeoExtent}
\alias{GeoExtent-class}
\alias{Rcpp_GeoExtent-class}

\alias{show,GeoExtent-method}

\title{Class "GeoExtent" }

\description{
Objects of class GeoExtent are used to define the spatial extent (extremes) of objects of the SpatRaster class. 
}

\section{Objects from the Class}{
You can use the \code{\link{ext}} function to create GeoExtent objects, or to extract them from SpatRaster objects.
}

\section{Slots}{
	 \describe{
    \item{\code{ptr}:}{pointer to the C++ class} 
  }
}

\section{Methods}{
  \describe{
 \item{show}{display values of a GeoExtent object }
 }
}


\examples{
e <- ext(-180, 180, -90, 90)
e
}

\keyword{classes}
\keyword{spatial}
