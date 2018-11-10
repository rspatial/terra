\name{writeRaster}

\alias{writeRaster,SpatRaster,character-method}

\alias{writeRaster}

\title{Write raster data to a file}

\description{
Write a SpatRaster object to a file. 
}

\usage{
\S4method{writeRaster}{SpatRaster,character}(x, filename, format="", datatype="FLT4S", overwrite=FALSE, ...)
}

\arguments{
\item{x}{Raster* object}
\item{filename}{Character. Output filename. Optional}
\item{format}{Character. Output file format. See \code{\link{writeRaster}}}
\item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
\item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{Additional arguments. None implemented}
}

\value{
This function is used for the side-effect of writing values to a file.
}


\keyword{ spatial }
\keyword{ methods }

