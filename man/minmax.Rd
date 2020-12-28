\name{range}

\alias{minmax}
\alias{setMinMax}
\alias{minmax,SpatRaster-method}
\alias{setMinMax,SpatRaster-method}

\title{Get or compute the minimum and maximum cell values}

\description{
The minimum and maximum value of a SpatRaster are returned or computed (from a file on disk if necessary) and stored in the object. 
}

\usage{
\S4method{minmax}{SpatRaster}(x)
\S4method{setMinMax}{SpatRaster}(x, force=FALSE)
}

\arguments{
  \item{x}{ SpatRaster }
  \item{force}{logical. If \code{TRUE} min and max values are recomputed even if already avaialbe }
}

\value{
setMinMax: nothing. Used for the side-effect of computing the minimum and maximum values of a SpatRaster

minmax: numeric matrix of minimum and maximum cell values by layer
}


\examples{ 
r <- rast(system.file("ex/elev.tif", package="terra"))
minmax(r)
}

\keyword{ spatial }
\keyword{ methods }

