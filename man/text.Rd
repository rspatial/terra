\name{text}

\docType{methods}

\alias{text}
\alias{text,SpatRaster-method}
\alias{text,SpatVector-method}

\title{Add labels to a map}

\description{
Plots labels, that is a textual (rather than color) representation of values, on top an existing plot (map). 
}

\usage{
\S4method{text}{SpatRaster}(x, labels, digits=0, halo=FALSE, hc="white", hw=0.1, jitter=0, ...)

\S4method{text}{SpatVector}(x, labels, halo=FALSE, inside=FALSE, hc="white", hw=0.1, jitter=0, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{labels}{character. Optional. Vector of labels with \code{length(x)} or a variable name from \code{names(x)}}
  \item{digits}{integer. How many digits should be used?}
  \item{halo}{logical. If \code{TRUE} a "halo" is printed around the text}
  \item{hc}{character. The halo color}
  \item{hw}{numeric. The halo width}
  \item{inside}{logical. Should the text always be placed inside one the sub-geometries?}
  \item{jitter}{numeric. Passed to \code{\link{jitter}} 'factor' argument. If >0, labels are randomly shifted from their centroid positions, possibly avoiding overlaps}
  \item{...}{additional arguments to pass to graphics function \code{\link[graphics]{text}} }
}

\seealso{
\code{\link[graphics]{text}, \link{plot}, \link{halo}}
}



\examples{ 
r <- rast(nrows=4, ncols=4)
values(r) <- 1:ncell(r)

plot(r)
text(r)

set.seed(123)
text(r, jitter = 2, col = "red", halo = TRUE)

plot(r)
text(r, halo=TRUE, hc="blue", col="white", hw=0.2)

plot(r, col=rainbow(16))
text(r, col=c("black", "white"), vfont=c("sans serif", "bold"), cex=2)
}

\keyword{methods}
\keyword{spatial}

