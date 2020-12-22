\name{fill}

\docType{methods}

\alias{fill}
\alias{fill,SpatVector-method}


\title{Remove holes from polygons}

\description{
Remove the holes in SpatVector polygons. If \code{inverse=TRUE} the holes are returned (as polygons).
} 

\usage{
\S4method{fill}{SpatVector}(x, inverse=FALSE, ...)
}

\arguments{
\item{x}{SpatVector}
\item{inverse}{logical. If \code{TRUE} the holes are returned as polygons}
\item{...}{additional arguments. None implemented }  
}

\value{SpatVector}

\examples{
x <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
hole <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))

z <- rbind(cbind(object=1, part=1, x, hole=0), 
		   cbind(object=1, part=1, hole, hole=1))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)
p <- vect(z, "polygons", atts=data.frame(id=1))
p

f <- fill(p)
g <- fill(p, inverse=TRUE)

plot(p, lwd=16, border="gray", col="light yellow")
polys(f, border="blue", lwd=3, density=4, col="orange")
polys(g, col="white", lwd=3)
}


\keyword{methods}
\keyword{spatial}


