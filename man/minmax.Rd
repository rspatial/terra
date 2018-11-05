\name{range}

\alias{minmax}
\alias{setMinMax}
\alias{minmax,SpatRaster-method}
\alias{setMinMax,SpatRaster-method}

\title{get or compute min and max values}

\description{
The minimum and maximum value of a SpatRaster are returned or computed (from a file on disk if necessary) and stored in the object. 
}

\usage{
\S4method{minmax}{SpatRaster}(x)
\S4method{setMinMax}{SpatRaster}(x)
}

\arguments{
  \item{x}{ SpatRaster }
}

\value{
setMinMax: nothing. Used for the side-effect of computing the minimum and maximum values of a SpatRaster

minmax: numeric matrix of minimum and maximum cell values by layer
}


\examples{ 
r <- rast(system.file("exdata/test.tif", package="terra"))
minmax(r)
}

\keyword{ spatial }
\keyword{ methods }

