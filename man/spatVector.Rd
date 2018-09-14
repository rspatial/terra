\name{spatVector}

\docType{methods}

\alias{spatPolygons}
\alias{spatPolygons,data.frame-method}

\title{Create SpatVector objects}

\description{
Create a new SpatPolygons object
} 

\usage{
\S4method{spatPolygons}{data.frame}(x, attr=NULL, crs="", ...)
}

\arguments{
\item{x}{data.frame}
\item{attr}{data.frame with the attributes. The number of rows must match the number of lines/polgyons}
\item{crs}{the coordinate reference system (PROJ4 notation)}
\item{...}{additional matrices and/or lists with matrices}  
}


\value{SpatPolygons object}


\examples{
x1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x4 <- rbind(c(41,-41.5), c(51,-35), c(62,-41), c(51,-50))

z <- rbind(cbind(object=1, part=1, x1), cbind(object=1, part=2, x2), cbind(object=2, part=1, x3), cbind(object=3, part=1,  x4))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)
z$hole <- z$part - 1

a <- spatPolysgons(z)
}

\keyword{methods}
\keyword{spatial}

