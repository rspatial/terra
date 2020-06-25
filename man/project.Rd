\name{project}

\alias{project}
\alias{project,SpatVector-method}
\alias{project,SpatRaster-method}

\title{Change the coordinate reference system}

\description{
Change the coordinate reference system ("project") of a SpatVector or SpatRaster. 
}


\usage{
\S4method{project}{SpatVector}(x, y, ...)

\S4method{project}{SpatRaster}(x, y, method="bilinear", mask=FALSE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster or character if \code{x} is a SpatRaster; character if \code{x} is a SpatVector. The character should be a PROJ.4 description of the output coordinate reference system (crs). If \code{x} is a SpatVector, you can also provide another object from which a crs can be extracted with \code{\link{crs}}}
  \item{method}{character. Method used for estimating the new cell values. One of: 
  
	near: nearest neighbour. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for continuous values.
	bilinear: bilinear interpolation. Default.
	cubic: cubic interpolation.
	cubicspline: cubic spline interpolation.
  }

  \item{mask}{logical. If \code{TRUE}, mask out areas outside the input extent (see example with Robinson projection)}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}   
}


\value{
SpatVector or SpatRaster
}

\seealso{ \code{\link{resample}}}


\examples{
# SpatRaster
a <- rast(ncol=40, nrow=40, xmin=-110, xmax=-90, ymin=40, ymax=60, 
          crs="+proj=longlat +datum=WGS84")
values(a) = 1:ncell(a)
newcrs="+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
b <- rast(ncol=94, nrow=124, xmin=-944881, xmax=935118, ymin=4664377, ymax=7144377, crs=newcrs)
w <- project(a, b)

if  (gdal_version() >= "3.0.0") {
# y can also be a crs
	ww <- project(a, newcrs)


# Global Robinson
	r <- rast()
	values(r) <- 1:ncell(r)
	crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
	x <- project(r, crs)
	# the corners should be empty, but are not
    plot(x)
	y <- project(r, crs, mask=TRUE)
    # now they are
	plot(y)
}




# SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"
p <- project(v, crs)
p
}

\keyword{spatial}

