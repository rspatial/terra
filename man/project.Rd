\name{project}

\alias{project}
\alias{project,SpatVector-method}

\title{Change the coordinate reference system of a SpatVector}

\description{
Change the coordinate reference system ("project") of a SpatVector
}


\usage{
\S4method{project}{SpatVector}(x, crs, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{crs}{character. Description of output coordinate reference system in PROJ.4 notation }
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
}

\keyword{spatial}

