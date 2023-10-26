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
\S4method{text}{SpatRaster}(x, labels, digits=0, halo=FALSE, ...)

\S4method{text}{SpatVector}(x, labels, halo=FALSE, inside=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{labels}{character. Optional. Vector of labels with \code{length(x)} or a variable name from \code{names(x)}}
  \item{digits}{integer. How many digits should be used?}
  \item{halo}{logical. If \code{TRUE} a "halo" is printed around the text. If \code{TRUE}, additional arguments \code{hc="white"} and \code{hw=0.1} can be modified to set the color and width of the halo}
  \item{inside}{logical. Should the text always be placed inside one the sub-geometries?}
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

plot(r)
text(r, halo=TRUE, hc="blue", col="white", hw=0.2)

plot(r, col=rainbow(16))
text(r, col=c("black", "white"), vfont=c("sans serif", "bold"), cex=2)
}

\keyword{methods}
\keyword{spatial}

