\name{extract}

\docType{methods}

\alias{extract}
\alias{extract,SpatRaster,SpatVector-method}
\alias{extract,SpatRaster,matrix-method}

\alias{[,SpatRaster,SpatVector,missing-method}
\alias{[,SpatRaster,missing,missing-method}
\alias{[,SpatRaster,numeric,missing-method}
\alias{[,SpatRaster,missing,numeric-method}
\alias{[,SpatRaster,numeric,numeric-method}


\title{Extract values from a SpatRaster}

\description{
Extract values from a SpatRaster for the locations of a SpatVector
} 

\usage{
\S4method{extract}{SpatRaster,SpatVector}(x, y, fun=NULL, ..., method="simple", drop=FALSE)
\S4method{extract}{SpatRaster,matrix}(x, y, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{y}{SpatVector (for points, lines, polygons), or 2-column matrix (x, y; for points)}
\item{fun}{character. function to summarize the data. Currently ignored}
\item{...}{additional arguments. None implemented}
\item{method}{character. method for extracting values with points. The default is "simple", the alternative is "bilinear"}
\item{drop}{boolean. If \code{TRUE}, the output is simplified}
}

\value{data.frame}

\seealso{\link{values}}

\examples{
r <- rast(ncol=5, nrow=5, xmin=0, xmax=5, ymin=0, ymax=5)
values(r) <- 1:25
xy <- rbind(c(0.5,0.5), c(2.5,2.5))
p <- vect(xy)

extract(r, xy)
extract(r, p)

r[1,]
r[5]
r[,5]

r[c(0:2, 99:101)]

f <- system.file("exdata/test.tif", package="terra")
r <- rast(f)

xy <- cbind(179000, 330000)
xy <- rbind(xy-100, xy, xy+1000)
extract(r, xy)

p <- vect(xy)
g <- geom(p)
g

extract(r, p)

x <- r + 10
extract(x, p)

i <- cellFromXY(r, xy)
x[i]
r[i]

y <- c(x,x*2,x*3)
y[i]

# extract with a polygon
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
z <- rast(v)
values(z) <- 1:100
e <- extract(z, v)
e[1:2]
rapply(e, mean)

x <- c(z, z*2, z/3)
e <- extract(x, v)
matrix(rapply(e, mean), ncol=nlyr(x), byrow=TRUE)
}

\keyword{methods}
\keyword{spatial}

