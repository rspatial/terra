\name{tags}

\alias{tags}
\alias{tags,SpatRaster-method}
\alias{set.tags}
\alias{set.tags,SpatRaster-method}

  
\title{Set or get tags}

\description{
You can set arbitrary metadata to a SpatRaster using "name=value" tags. 
}

\usage{
\S4method{set.tags}{SpatRaster}(x, name, value)

\S4method{tags}{SpatRaster}(x, name=NULL)
}

\arguments{
\item{x}{SpatRaster}

\item{name}{character}

\item{value}{character}
}

\value{
SpatRaster (\code{set.tags}), matrix (\code{tags}) or character (\code{tags} and \code{name} is not NULL)
}

\examples{
r <- rast()
r <- set.tags(r, c("one", "two", "three"), c("ABC", "123", "hello"))
tags(r)

r <- set.tags(r, "another_tag", "another_value")
tags(r)

tags(r, "two")

# remove one tag
r <- set.tags(r, c("one"), c(""))
tags(r)

# remove all tags
r <- set.tags(r, NULL)
tags(r)
}

\keyword{spatial}

