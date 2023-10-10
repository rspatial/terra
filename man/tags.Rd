\name{tags}

\alias{tags}
\alias{tags,SpatRaster-method}
\alias{set.tags}
\alias{set.tags,SpatRaster-method}

  
\title{Set or get tags to a SpatRaster}

\description{
You can set arbitrary tags to a SpatRaster. Each tag has a name and a value.
}

\usage{
\S4method{set.tags}{SpatRaster}(x, name, value)

\S4method{tags}{SpatRaster}(x)
}

\arguments{
\item{x}{SpatRaster}

\item{name}{character}

\item{value}{character}
}

\value{
SpatRaster (\code{set.tags}) or matrix (\code{tags})
}

\examples{
r <- rast()
r <- set.tags(r, c("one", "two", "three"), c("ABC", "123", "hello"))
tags(r)

r <- set.tags(r, "another_tag", "another_value")
tags(r)

# remove one tag
r <- set.tags(r, c("one"), c(""))
tags(r)

# remove all tags
r <- set.tags(r, NULL)
tags(r)
}

\keyword{spatial}

