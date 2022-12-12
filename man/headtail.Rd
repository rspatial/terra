\name{headtail}

\docType{methods}

\alias{head}
\alias{head,SpatRaster-method}
\alias{head,SpatVector-method}
\alias{tail}
\alias{tail,SpatRaster-method}
\alias{tail,SpatVector-method}

\title{head and tail of a SpatRaster or SpatVector}

\description{
Show the head (first values) or tail (last values) of a SpatRaster or of the attributes of a SpatVector. 
}

\usage{
head(x, ...) 
tail(x, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{...}{additional arguments passed on to other methods}
}

\value{
matrix (SpatRaster) or data.frame (SpatVector)
}

\seealso{
\code{\link{show}}, \code{\link{geom}}
}



\examples{
r <- rast(nrows=25, ncols=25)
values(r) <- 1:ncell(r)
head(r)
tail(r)
}

\keyword{methods}
\keyword{spatial}

