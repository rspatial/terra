\name{tags}

\alias{tags}
\alias{tags,SpatRaster-method}
\alias{tags<-}
\alias{tags<-,SpatRaster-method}

  
\title{Set or get tags}

\description{
You can set arbitrary metadata to a SpatRaster using "name=value" tags. 
}

\usage{
\S4method{tags}{SpatRaster}(x)<-value

\S4method{tags}{SpatRaster}(x, name=NULL)
}

\arguments{
\item{x}{SpatRaster}

\item{name}{character}

\item{value}{character of "name=value" or two-column matrix}
}

\value{
SpatRaster (\code{tags<-}), or named character (\code{tags})
}

\examples{
r <- rast(ncol=5, nrow=5)
m <- cbind(c("one", "two", "three"), c("ABC", "123", "hello"))
tags(r) <- m
tags(r)

tags(r) <- c("another_tag=another_value", "one more=this value")
tags(r)

tags(r, "two")

# remove a tag
tags(r) <- cbind("one", "")
tags(r) <- "two="
tags(r)

# remove all tags
tags(r) <- NULL
tags(r)
}

\keyword{spatial}

