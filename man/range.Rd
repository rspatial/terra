\name{range}

\alias{range}
\alias{setRange}
\alias{range,SpatRaster-method}
\alias{setRange,SpatRaster-method}

\title{get or compute min and max values}

\description{
The minimum and maximum value of a SpatRaster are computed (from a file on disk if necessary) and stored in the object. 
}

\usage{
\S4method{setRange}{SpatRaster}(x)
\S4method{range}{SpatRaster}(x, ..., na.rm=TRUE)
}

\arguments{
  \item{x}{ SpatRaster object }
  \item{...}{ additional arguments. None implemented }
  \item{na.rm}{ ignored }
}

\value{
getRange: nothing. Used for the side-effect of setting the range on a SpatRaster object

range: numeric matrix of min and max values
}


\examples{ 
r <- georst(system.file("external/test.grd", package="terra"))
range(r)
}

\keyword{ spatial }
\keyword{ methods }

