\name{bestMatch}

\alias{bestMatch}
\alias{bestMatch,SpatRaster,SpatVector-method}
\alias{bestMatch,SpatRaster,data.frame-method}

\title{bestMatch}

\description{
Determine for each grid cell which reference it is most similar to. A reference consists of a SpatVector with reference locations, or a data.frame in which each column matches a layer name in the SpatRaster.

Similarity is computed with the sum of squared differences between the cell and the reference. It may be important to first scale the input. 
}

\usage{
\S4method{bestMatch}{SpatRaster,SpatVector}(x, y, labels="", filename="", ...)

\S4method{bestMatch}{SpatRaster,data.frame}(x, y, labels="", filename="", ...)
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
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)

# locations of interest 
pts <- vect(cbind(c(25.25, 34.324, 43.003), c(54.577, 46.489, 30.905)))
pts$code <- LETTERS[1:3]

plot(r)
points(pts, pch=20, cex=2, col="red")
text(pts, "code", pos=4, halo=TRUE)

x <- scale(r)

s1 <- bestMatch(x, pts, labels=pts$code)
plot(s1)

# same result
e <- extract(x, pts, ID=FALSE)
s2 <- bestMatch(x, e, labels=c("Ap", "Nt", "Ms"))
}

\keyword{spatial}

