\name{hist}

\docType{methods}

\alias{hist}
\alias{hist,SpatRaster-method}

\title{Histogram}

\description{
Create a histogram of the values of a SpatRaster. For large datasets a sample of \code{maxcell} is used.
}

\usage{
\S4method{hist}{SpatRaster}(x, layer, maxcell=1000000, plot=TRUE, maxnl=16, main, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer or character to indicate layer numbers (or names). If missing, all layers up to \code{maxnl} are used} 
  \item{maxcell}{integer. To regularly sample very large objects}
  \item{plot}{logical. Plot the histogram or only return the histogram values}
  \item{maxnl}{positive integer. The maximum number of layers to use. Ignored if \code{layer} is not missing}
  \item{main}{character. Main title(s) for the plot. Default is the value of \code{\link{names}}}
  \item{...}{additional arguments. See \code{\link[graphics]{hist}}}  
}


\value{
This function is principally used for plotting a histogram, but it also returns an object of class "histogram" (invisibly if \code{plot=TRUE}).
}

\seealso{ \code{\link{pairs}, \link{boxplot}} }

\examples{ 
r1 <- r2 <- rast(nrows=50, ncols=50)
values(r1) <- runif(ncell(r1))
values(r2) <- runif(ncell(r1))
rs <- r1 + r2
rp <- r1 * r2

opar <- par(no.readonly =TRUE)
par(mfrow=c(2,2))
plot(rs, main='sum')
plot(rp, main='product')
hist(rs)
a <- hist(rp)
a
x <- c(rs, rp, sqrt(rs))
hist(x)
par(opar)
}

\keyword{methods}
\keyword{spatial}

