\name{as.character}

\docType{methods}

\alias{as.character,SpatExtent-method}
\alias{as.character,SpatRaster-method}


  
\title{Create a text representation of (the skeleton of) an object}

\description{
Create a text representation of (the skeleton of) an object
}

\usage{
\S4method{as.character}{SpatExtent}(x)

\S4method{as.character}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
}
 
\value{
character
}
  

\examples{
r <- rast()
ext(r)
ext(c(0, 20, 0, 20))

}

\keyword{spatial}

