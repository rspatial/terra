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
\S4method{text}{SpatRaster}(x, labels, digits=0, halo=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...)

\S4method{text}{SpatVector}(x, labels, halo=FALSE, inside=FALSE, fit=FALSE, hc="white", hw=0.1, jitter=0, xpd=TRUE, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{labels}{character. Optional. Vector of labels with \code{length(x)} or a variable name from \code{names(x)}}
  \item{digits}{integer. How many digits should be used?}
  \item{halo}{logical. If \code{TRUE} a "halo" is printed around the text}
  \item{hc}{character. The halo color}
  \item{hw}{numeric. The halo width}
  \item{inside}{logical. if \(x) is a SpatVector of polygons, should the (x,y) coordinates of the label be inside the geometry}
  \item{fit}{logical. if \(x) is a SpatVector of polygons, should an attempt be made assure that the labels are entirely inside their corresponding polygon}
  \item{jitter}{numeric. The amount of random noise used to adjust label positions, possibly avoiding overlaps. See argument 'factor' in \code{\link{jitter}}}
  \item{xpd}{logical or \code{NA}. If \code{FALSE}, the text is clipped to the plot region, if \code{TRUE}, the text is clipped to the figure region, and if \code{NA}, the text is clipped to the device region}
  \item{...}{additional arguments to pass to graphics function \code{\link[graphics]{text}} }
}

\seealso{
\code{\link[graphics]{text}, \link{plot}, \link{halo}}
}

\value{
matrix with the label locations (invisibly)
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
xy <- text(r, col=c("black", "white"), vfont=c("sans serif", "bold"), cex=2)
head(xy)
}

\keyword{methods}
\keyword{spatial}

