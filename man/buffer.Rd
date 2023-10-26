\name{buffer}

\alias{buffer}
\alias{buffer,SpatRaster-method}
\alias{buffer,SpatVector-method}


\title{Create a buffer around vector geometries or raster patches}

\description{
Calculate a buffer around all cells that are not \code{NA} in a SpatRaster, or around the geometries of a SpatVector.

SpatRaster cells inside the buffer distance get a value of 1.
 
Note that the distance unit of the buffer \code{width} parameter is meters if the CRS is (\code{+proj=longlat}), and in map units (typically also meters) if not. 
}

\usage{
\S4method{buffer}{SpatRaster}(x, width, background=0, filename="", ...)

\S4method{buffer}{SpatVector}(x, width, quadsegs=10, capstyle="round", 
		joinstyle="round", mitrelimit=NA, singlesided=FALSE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{width}{numeric. Unit is meter if \code{x} has a longitude/latitude CRS, or in the units of the coordinate reference system in other cases (typically also meter). The value should be > 0 if \code{x} is a SpatRaster. If \code{x} is a SpatVector, this argument is vectorized, meaning that you can provide a different value for each geometry in \code{x}; and you can also use the name of a variable in \code{x} that has the widths}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{background}{numeric. value to assign to cells outside the buffer. If this value is zero or FALSE, a boolean SpatRaster is returned}  
  \item{quadsegs}{positive integer. Number of line segments to use to draw a quart circle}  
  \item{capstyle}{character. One of "round", "square" or "flat". Ignored if \code{is.lonlat(x)}}
  \item{joinstyle}{character. One of "round", "mitre" or "bevel". Ignored if \code{is.lonlat(x)}}  
  \item{mitrelimit}{numeric. Place an upper bound on a mitre join to avoid it from extending very far from acute angles in the input geometry. Ignored if \code{is.lonlat(x)}}
  \item{singlesided}{logical. If \code{TRUE} a buffer is constructed on only one side of each input line.  Ignored if \code{is.lonlat(x)}}  
}


\value{
Same as \code{x} 
}


\seealso{
\code{\link{distance}}, \code{\link{elongate}}
}


\examples{
r <- rast(ncols=36, nrows=18)
r[500] <- 1
b <- buffer(r, width=5000000) 
plot(b)

v <- vect(rbind(c(10,10), c(0,60)), crs="+proj=merc")
b <- buffer(v, 20)
plot(b)
points(v)

crs(v) <- "+proj=longlat" 
b <- buffer(v, 1500000)
plot(b)
points(v)
}

\keyword{spatial}

