\name{extremes}

\alias{minmax}
\alias{hasMinMax}
\alias{setMinMax}
\alias{minmax,SpatRaster-method}
\alias{hasMinMax,SpatRaster-method}
\alias{setMinMax,SpatRaster-method}

\title{Get or compute the minimum and maximum cell values}

\description{
The minimum and maximum value of a SpatRaster are returned or computed (from a file on disk if necessary) and stored in the object. 
}

\usage{
\S4method{minmax}{SpatRaster}(x, compute=FALSE)
\S4method{hasMinMax}{SpatRaster}(x)
\S4method{setMinMax}{SpatRaster}(x, force=FALSE)
}

\arguments{
  \item{x}{ SpatRaster }
  \item{compute}{logical. If \code{TRUE} min and max values are computed if they are not available}
  \item{force}{logical. If \code{TRUE} min and max values are recomputed even if already available}
}

\value{
minmax: numeric matrix of minimum and maximum cell values by layer

hasMinMax: logical indicating whether the min and max values are available.

setMinMax: nothing. Used for the side-effect of computing the minimum and maximum values of a SpatRaster
}


\examples{ 
r <- rast(system.file("ex/elev.tif", package="terra"))
minmax(r)
}

\keyword{ spatial }
\keyword{ methods }

