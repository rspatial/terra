\name{dots}

\docType{methods}

\alias{dots}
\alias{dots,SpatVector-method}


\title{Make a dot-density map}

\description{
Create the dots for a dot-density map and add these to the current map. Dot-density maps are made to display count data. For example of population counts, where each dot represents n persons. The dots are returned as a \code{SpatVector}. It there is an active graphics device, the dots are added to it with \code{\link{points}}.
}

\usage{
\S4method{dots}{SpatVector}(x, field, size,  ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{field}{character of numeric indicating field name. Or numeric vector of the same length as \code{x}}
  \item{size}{positive number indicating the number of cases associated with each dot}
  \item{...}{graphical arguments passed to \code{points}}
}


\value{
SpatVector (invisibly)
}
 
\seealso{
 \code{\link{plot}}, \code{\link{cartogram}}, \code{\link{points}}
}

\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$population <- 1000*(1:12)^2
plot(v, lwd=3, col="light gray", border="white")
d <- dots(v, "population", 1000, col="red", cex=.75)
lines(v)
d
}

\keyword{methods}
\keyword{spatial}

