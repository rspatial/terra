\name{writeVector}

\alias{writeVector,SpatVector,character-method}

\alias{writeVector}

\title{Write vector data to a file}

\description{
Write a SpatVector object to a file. 
}

\usage{
\S4method{writeVector}{SpatVector,character}(x, filename, overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files. See Details }
  \item{...}{additional arguments. None implemented}
}

\value{
SpatVector (invisibly). This function is used for the side-effect of writing values to a file.
}


\keyword{ spatial }
\keyword{ methods }

