\docType{methods}
\name{pitfinder2}

\alias{pitfinder2}
\alias{pitfinder2,SpatRaster-method}


\title{Pit Finder in a Flow Dir Raster map for Waterhed Extraction}

\description{
It extractions a watershed mask from a flow direction input raster point a pour-point cell.
}

\usage{
\S4method{pitfinder2}{SpatRaster}(p,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 for the pits  and value 0 elsewhere. 
}

\details{
Lorem Ipsum. You can insert information about tha algorithm
}


\author{

Emanuele Cordano
}

\examples{


## TO DO


}


\references{
Lorem Ipsum
}



\keyword{spatial}
