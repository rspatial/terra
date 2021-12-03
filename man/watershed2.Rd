\name{watershed2}

\alias{watershed2}
\alias{watershed2,SpatRaster-method}

\title{Waterhed Extraction}

\description{
It extractions a watershed mask from a flow direction input raster point a pour-point cell.
}

\usage{
\S4method{watershed2}{SpatRaster}(p,pp_offset,filename, filename="",...)  
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{pp_offset}{pour point offset.}
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 in the extracted watershed mask and value 0 elsewhere. 
}

\details{
Lorem Ipsum. You can insert information about tha algorithm
}


\author{

Ezio Crestaz, Emanuele Cordano
}

\references{
Lorem Ipsum
}



\keyword{spatial}
