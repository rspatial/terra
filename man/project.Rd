\name{project}

\alias{project}
\alias{project,SpatVector-method}

\title{Change the coordinate reference system}

\description{
Change the coordinate reference system ("project") of a SpatVector or SpatRaster
}


\usage{
\S4method{project}{SpatVector}(x, crs, ...)
\S4method{project}{SpatRaster}(x, crs, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{crs}{character. Description of output coordinate reference system in PROJ.4 notation }
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}   
  \item{...}{Additional arguments. None implemented}  
}


\value{
SpatVector 
}

\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"
p <- project(v, crs)
p

a <- rast(ncol=40, nrow=40, xmin=-110, xmax=-90, ymin=40, ymax=60, crs="+proj=longlat +datum=WGS84")
values(a) = 1:ncell(a)
crs="+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
b <- project(a, crs)
}

\keyword{spatial}

