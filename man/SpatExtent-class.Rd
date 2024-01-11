\name{SpatExtent-class}

\docType{class}

\alias{SpatExtent}
\alias{SpatExtent-class}
\alias{Rcpp_SpatExtent-class}

\alias{show,SpatExtent-method}

\title{Class "SpatExtent" }

\description{
Objects of class SpatExtent are used to define the spatial extent (extremes) of objects of the SpatRaster class. 
}

\section{Objects from the Class}{
You can use the \code{\link{ext}} function to create SpatExtent objects, or to extract them from a SpatRaster, SpatVector or related objects.
}


\section{Methods}{
  \describe{
 \item{show}{display values of a SpatExtent object }
 }
}


\examples{
e <- ext(-180, 180, -90, 90)
e
}

\keyword{classes}
\keyword{spatial}
