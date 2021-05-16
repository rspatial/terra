\name{mem}

\alias{mem}

\title{Memory available and needed}

\description{
\code{mem_info} prints the amount of RAM that is required and available to process a SpatRaster. 

\code{free_RAM} returns the amount of RAM that is available
}

\usage{
mem_info(x, n=1, print=TRUE)

free_RAM()
}

\arguments{
  \item{x}{SpatRastaer}
  \item{n}{positive integer. The number of copies of \code{x} that are needed}
  \itme{print}{logical. Print information to console}
}


\value{
free_RAM returns the amount of available RAM in kilobytes
}

\examples{
mem_info(rast())

free_RAM()
}

\keyword{spatial}

