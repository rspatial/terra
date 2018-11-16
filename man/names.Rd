\name{names}

\alias{names}
\alias{names<-}
\alias{names,SpatRaster-method}
\alias{names<-,SpatRaster-method}
\alias{names,SpatVector-method}
  
\title{Names of Spat* objects}

\description{
Get or set the names of the layers of a SpatRaster.
}

\usage{
\S4method{names}{SpatRaster}(x)
\S4method{names}{SpatRaster}(x)<-value
}

\arguments{
  \item{x}{For names: SpatRaster. For validNames: character}
  \item{value}{character (vector)}
}

\value{
Character 
}


\examples{
s <- rast(ncols=5, nrows=5, nlyr=3)
nlyr(s)
names(s)
names(s) <- c('a', 'b', 'c')
# space is not valid
names(s)[2] <- 'hello world'
names(s)

# two invalid names
names(s) <- c('a', '  a  ', '3')
s

}

\keyword{spatial}
