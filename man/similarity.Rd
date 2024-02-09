\name{similarity}

\alias{similarity}
\alias{similarity,SpatRaster,SpatVector-method}
\alias{similarity,SpatRaster,data.frame-method}

\title{Similarity}

\description{
Determine for each grid cell the location/record it is most similar to. 
}

\usage{
\S4method{similarity}{SpatRaster,SpatRaster}(x, y, labels="", filename="", ...)

\S4method{similarity}{SpatRaster,data.frame}(x, y, labels="", filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or data.frame}
  \item{labels}{character. labels that correspond to each class (row in \code{y}}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}


\examples{
library(terra)
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)

# locations of interest 
pts <- vect(cbind(c(25.25, 34.324, 43.003), c(54.577, 46.489, 30.905)))
pts$code <- LETTERS[1:3]

plot(r)
points(pts, pch=20, cex=2, col="red")
text(pts, "code", pos=4, halo=TRUE)

x <- scale(r)

s1 <- similarity(x, pts, labels=pts$code)
plot(s1)

# same result
e <- extract(x, pts, ID=FALSE)
s2 <- similarity(x, e, labels=c("Ap", "Nt", "Ms"))
}

\keyword{spatial}

