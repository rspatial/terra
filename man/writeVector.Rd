\name{writeVector}

\alias{writeVector,SpatVector,character-method}

\alias{writeVector}

\title{Write SpatVector data to a file}

\description{
Write a SpatVector to a file. You can choose one of many file formats.  
}

\usage{
\S4method{writeVector}{SpatVector,character}(x, filename, filetype=NULL, layer=NULL, overwrite=FALSE, options="ENCODING=UTF-8")
}

\arguments{
  \item{x}{SpatVector}
  \item{filename}{character. Output filename}
  \item{filetype}{character. A file format associated with a GDAL "driver" such as "ESRI Shapefile". See \code{ gdal(drivers=TRUE)}. If \code{NULL} it is attempted to guess the filetype from the filename extension}
  \item{layer}{character. Output layer name. If \code{NULL} the filename is used}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{options}{character. Format specific GDAL options such as "ENCODING=UTF-8". Use NULL or "" to not use any options}
}


\examples{
v <- vect(cbind(1:5,1:5))
crs(v) <- "+proj=longlat +datum=WGS84"
v$id <- 1:length(v)
v$name <- letters[1:length(v)]
tmpf1 <- tempfile()
writeVector(v, tmpf1)
x <- vect(tmpf1)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
tmpf2 <- tempfile()
writeVector(v, tmpf2)
y <- vect(tmpf2)
}



\keyword{ spatial }
\keyword{ methods }

