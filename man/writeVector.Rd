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


\keyword{ spatial }
\keyword{ methods }

