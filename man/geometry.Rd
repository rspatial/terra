\name{geom}

\docType{methods}

\alias{geom}
\alias{geom,SpatLayer-method}

\title{geom}

\description{
Get the geom of a SpatLayer
}

\usage{
\S4method{geom}{SpatLayer}(x, ...)

}

\arguments{
  \item{x}{SpatLayer}
  \item{...}{Additional arguments. None implemented}  
}


\value{
matrix
}


\examples{
x1 <- rbind(c(-175,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x4 <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=0, part=0, x1), cbind(object=1, part=0, x2), 
           cbind(object=2, part=0, x3), cbind(object=2, part=1,  x4))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)
z$hole <- 0
z$hole[z$object==2 & z$part==1] <- 1

p <- vect(z)
geom(p)

f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
g <- geom(v)
head(g)
}

\keyword{spatial}
