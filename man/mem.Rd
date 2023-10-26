\name{mem}

\alias{mem_info}
\alias{free_RAM}


\title{Memory available and needed}

\description{
\code{mem_info} prints the amount of RAM that is required and available to process a SpatRaster. 

\code{free_RAM} returns the amount of RAM that is available
}

\usage{
mem_info(x, n=1)

free_RAM()
}

\arguments{
  \item{x}{SpatRaster}
  \item{n}{positive integer. The number of copies of \code{x} that are needed}
}


\value{
free_RAM returns the amount of available RAM in kilobytes
}

\examples{
mem_info(rast())

free_RAM()
}

\keyword{spatial}

