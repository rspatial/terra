\name{vector_layers}

\alias{vector_layers}

\title{List or remove layers from a vector file}

\description{
List or remove layers from a vector file that supports layers such as GPGK
}

\usage{
vector_layers(filename, delete="", return_error=FALSE)
}

\arguments{
  \item{filename}{character. filename}
  \item{delete}{character. layers to be deleted (ignored if the value is \code{""}}
  \item{return_error}{logical. If \code{TRUE}, an error occurs if some layers cannot be deleted. Otherwise a warning is given}
}



\keyword{ spatial }
\keyword{ methods }

