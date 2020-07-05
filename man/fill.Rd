\name{fill}

\docType{methods}

\alias{fill}
\alias{fill,SpatVector-method}


\title{Remove holes from polygons}

\description{
Remove the holes in SpatVector polygons
} 

\usage{
\S4method{fill}{SpatVector}(x, ...)
}

\arguments{
\item{x}{SpatVector}
\item{...}{additional arguments. None implemented }  
}

\value{character}

\examples{
x <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
hole <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))

z <- rbind(cbind(object=1, part=1, x, hole=0), 
		   cbind(object=1, part=1, hole, hole=1))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)
p <- vect(z, "polygons")
p

f <- fill(p)

plot(f, col="light blue", lwd=8, border="blue")
lines(p, lwd=2, col="red")
}


\keyword{methods}
\keyword{spatial}


