\name{scoff}

\alias{scoff}
\alias{scoff<-}
\alias{scoff,SpatRaster-method}
\alias{scoff<-,SpatRaster-method}


\title{Scale (gain) and offset}


\description{
These functions can be used to get or set the scale (gain) and offset parameters used to transform values when reading raster data from a file. The parameters are applied to the raw values using the formula below: 

\code{value <- value * scale + offset}

The default value for scale is 1 and for offset is 0. 'scale' is sometimes referred to as 'gain'.

Note that setting the scale and/or offset are intended to be used with values that are stored in a file. When values are memory, assigning scale or offset values will lead to the immediate computation of new values; in such cases it would be clearer to use \code{\link[terra]{Arith-methods}}.
}

\usage{
\S4method{scoff}{SpatRaster}(x)

\S4method{scoff}{SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{two-column matrix with scale (first column) and offset (second column) for each layer. Or \code{NULL} to remove all scale and offset values}
}

\value{
matrix or changed SpatRaster
}



\examples{ 
r <- rast(system.file("ex/elev.tif", package="terra"))
minmax(r)
scoff(r)
r[4603]

scoff(r) <- cbind(10, 5)
minmax(r)
scoff(r)
r[4603]
}

\keyword{ spatial }
\keyword{ methods }
