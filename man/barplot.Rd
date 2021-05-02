\name{barplot}

\docType{methods}

\alias{barplot}
\alias{barplot,SpatRaster-method}

\title{Bar plot of a SpatRaster}

\description{Create a barplot of the values of a the first layer of a SpatRaster. For large datasets a regular sample with a size of approximately \code{maxcells} is used.}

\usage{
\S4method{barplot}{SpatRaster}(height, maxcell=1000000, digits=0, breaks=NULL, col, ...)
}

\arguments{
  \item{height}{SpatRaster}
  \item{maxcell}{integer. To regularly subsample very large datasets}
  \item{digits}{integer used to determine how to \code{\link{round}} the values before tabulating. Set to \code{NULL} or to a large number if you do not want any rounding }
  \item{breaks}{breaks used to group the data as in \code{\link[base]{cut}}}
  \item{col}{a color generating function such as \code{\link{rainbow}} (the default), or a vector of colors}
  \item{...}{additional arguments for plotting as in \code{\link[graphics]{barplot}}}
}


\seealso{ \code{\link{hist}, \link{boxplot}} }

\value{
A numeric vector (or matrix, when \code{beside = TRUE}) of the coordinates of the bar midpoints, useful for adding to the graph. See \code{\link[graphics]{barplot}}
}


\examples{ 
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
barplot(r, digits=-1, las=2, ylab="Frequency")

op <- par(no.readonly = TRUE)
par(mai = c(1, 2, .5, .5))
barplot(r, breaks=10, col=c("red", "blue"), horiz=TRUE, digits=NULL, las=1)
par(op)
}

\keyword{methods}
\keyword{spatial}

