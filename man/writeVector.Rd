\name{writeVector}

\alias{writeVector,SpatVector,character-method}

\alias{writeVector}

\title{Write vector data to a file}

\description{
Write a SpatVector object to a file. 
}

\usage{
\S4method{writeVector}{SpatVector,character}(x, filename, overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments. None implemented}
}

\value{
SpatVector (invisibly). This function is used for the side-effect of writing values to a file.
}

\examples{
v <- vect(cbind(1:5,1:5))
crs(v) <- "+proj=lonlat +datum=WGS84"
v$id <- 1:length(v)
v$name <- letters[1:length(v)]
f <- tempfile()
writeVector(v, f)
x <- vect(f)
}



\keyword{ spatial }
\keyword{ methods }

