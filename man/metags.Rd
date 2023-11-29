\name{metags}

\alias{metags}
\alias{metags,SpatRaster-method}
\alias{metags<-}
\alias{metags<-,SpatRaster-method}

  
\title{Set or get metadata}

\description{
You can set arbitrary metadata to a SpatRaster using "name=value" tags. 
}

\usage{
\S4method{metags}{SpatRaster}(x, layer=NULL)<-value

\S4method{metags}{SpatRaster}(x, layer=NULL, name=NULL)
}

\arguments{
\item{x}{SpatRaster}
\item{layer}{NULL or positive integer. If the value is NULL, the tags assigned or returned are for the SpatRaster. Otherwise for the layer number(s)}

\item{name}{character}

\item{value}{character of "name=value" or two-column matrix}
}

\value{
SpatRaster (\code{metags<-}), or named character (\code{metags})
}

\examples{
r <- rast(ncol=5, nrow=5)
m <- cbind(c("one", "two", "three"), c("ABC", "123", "hello"))
metags(r) <- m
metags(r)

metags(r) <- c("another_tag=another_value", "one more=this value")
metags(r)

metags(r) <- c(another_tag="44", `one more`="that value")
metags(r)


metags(r, name="two")

# remove a tag
metags(r) <- cbind("one", "")
metags(r) <- "two="
metags(r)

# remove all metags
metags(r) <- NULL
metags(r)
}

\keyword{spatial}

