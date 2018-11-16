\name{head and tail}

\docType{methods}

\alias{head}
\alias{head,SpatRaster-method}
\alias{head,SpatVector-method}
\alias{tail}
\alias{tail,SpatRaster-method}
\alias{tail,SpatVector-method}

\title{Show the head or tail of a Spat* object}

\description{
Show the head (first values) or tail (last values) of a SpatRaster or of the attributes of a SpatVector. 
}

\usage{
head(x, ...) 
tail(x, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{...}{Additional arguments passed on to other methods}
}

\value{
matrix (SpatRaster) or data.frame (SpatVector)
}

\seealso{
\code{\link{show}}, \code{\link{geom}}
}



\examples{
r <- rast(nrow=25, ncol=25)
values(r) <- 1:ncell(r)
head(r)
tail(r)
}

\keyword{methods}
\keyword{spatial}

