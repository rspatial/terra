\name{boxplot}

\docType{methods}

\alias{boxplot}
\alias{boxplot,SpatRaster-method}

\title{
Box plot of SpatRaster data
}

\description{
Box plot of layers in a SpatRaster 
}

\usage{
\S4method{boxplot}{SpatRaster}(x, y=NULL, maxcell=100000, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{y}{NULL or a SpatRaster. If \code{x} is a SpatRaster it used to group the values of \code{x} by "zone"}
\item{maxcell}{Integer. Number of cells to sample from datasets}
\item{...}{additional arguments passed to \code{graphics::\link[graphics]{boxplot}}}
}

\value{
boxplot returns a list (invisibly) that can be used with \code{\link{bxp}}
}

\seealso{
\code{\link{pairs}, \link{hist}}
}

\examples{
r1 <- r2 <- r3 <- rast(ncols=10, nrows=10)
set.seed(409)
values(r1) <- rnorm(ncell(r1), 100, 40)
values(r2) <- rnorm(ncell(r1), 80, 10)
values(r3) <- rnorm(ncell(r1), 120, 30)
s <- c(r1, r2, r3)
names(s) <- c("Apple", "Pear", "Cherry")

boxplot(s, notch=TRUE, col=c("red", "blue", "orange"), main="Box plot", ylab="random", las=1)

op <- par(no.readonly = TRUE)
par(mar=c(4,6,2,2))
boxplot(s, horizontal=TRUE, col="lightskyblue", axes=FALSE)
axis(1)
axis(2, at=0:3, labels=c("", names(s)), las=1, cex.axis=.9, lty=0)
par(op)

## boxplot with 2 layers
v <- vect(system.file("ex/lux.shp", package="terra"))
r <- rast(system.file("ex/elev.tif", package="terra"))
y <- rasterize(v, r, "NAME_2")
b <- boxplot(r, y)
bxp(b)
}


\keyword{spatial}

