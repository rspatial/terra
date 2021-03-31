\name{merge}

\docType{methods}

\alias{merge}
\alias{merge,SpatRaster,SpatRaster-method}
\alias{merge,SpatExtent,SpatExtent-method}
\alias{merge,SpatVector,data.frame-method}

\title{ 
Merge SpatRaster or SpatExtent objects, or a SpatVector with a data.frame
}

\description{
Merge SpatRasters to form a new SpatRaster object with a larger spatial extent. If objects overlap, the values get priority in the same order as the arguments. See \code{\link{classify}} to merge a \code{SpatRaster} and a \code{data.frame}. You can also merge SpatExtent objects.

There is a also a method for merging SpatVector with a data.frame; that is, to join the data.frame to the attribute table of the SpatVector..
}

\usage{
\S4method{merge}{SpatRaster,SpatRaster}(x, y, ..., filename="", overwrite=FALSE, wopt=list())

\S4method{merge}{SpatExtent,SpatExtent}(x, y, ...)

\S4method{merge}{SpatVector,data.frame}(x, y, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatExtent}
  \item{y}{object of same class as \code{x}}
  \item{...}{if \code{x} is a SpatRaster or SpatVector: additional objects of the same class as \code{x}. If \code{x} is a SpatVector, the same arguments as in \code{\link[base]{merge}}}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 

\details{
The SpatRaster objects must have the same origin and spatial resolution. In areas where the SpatRaster objects overlap, the values of the SpatRaster that is last in the sequence of arguments will be retained. 
}

\value{
SpatRaster or SpatExtent
}

\seealso{
\code{\link{mosaic}}
}

\note{
You can use \code{merge} with \code{do.call} to merge a list of SpatRasters (see example). But note that if the list is named, these names are used by \code{merge}. So if all elements are named, there should be one element with a SpatRaster called \code{x} and another one called \code{y}. For example with \code{names(x)[1:2] <- c("x"m "y")}. You can also removed the names of the the first two elements (assuming these are SpatRasters) with \code{names(x)[1:2] <- ""}.
}


\examples{
x <- rast(xmin=-110, xmax=-50, ymin=40, ymax=70, ncols=60, nrows=30)
y <- rast(xmin=-80, xmax=-20, ymax=60, ymin=30)
res(y) <- res(x)
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)
mr <- merge(x, y)
plot(mr)
mr <- merge(y, x)

# if you have many SpatRaster objects in a list
# you can use do.call:
s <- list(x, y)
# add arguments such as filename
s$filename <- ""
m <- do.call(merge, s)

##
# SpatVector with data.frame
f <- system.file("ex/lux.shp", package="terra")
p <- vect(f)
dfr <- data.frame(District=p$NAME_1, Canton=p$NAME_2, Value=round(runif(length(p), 100, 1000)))
dfr <- dfr[1:5, ]
pm <- merge(p, dfr, all.x=TRUE, by.x=c('NAME_1', 'NAME_2'), by.y=c('District', 'Canton'))
pm
values(pm)
}

\keyword{methods}
\keyword{spatial}

