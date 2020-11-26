\name{cells}

\docType{methods}

\alias{cells}
\alias{cells,SpatRaster,missing-method}
\alias{cells,SpatRaster,vector-method}
\alias{cells,SpatRaster,SpatVector-method}
\alias{cells,SpatRaster,SpatExtent-method}


\title{Get cell numbers}

\description{ 
Get the cell numbers covered by a SpatVector or SpatExtent. Or that match values in a vector; or all non \code{NA} values.
}

\usage{
\S4method{cells}{SpatRaster,missing}(x, y, ...)

\S4method{cells}{SpatRaster,numeric}(x, y, ...)

\S4method{cells}{SpatRaster,SpatVector}(x, y, touches=is.lines(y), method="simple", ...)

\S4method{cells}{SpatRaster,SpatExtent}(x, y, ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector, SpatExtent, 2-column matrix representing points, numeric representing values to match, or missing}
  \item{touches}{logical. If \code{TRUE}, values for all cells touched by lines or polygons are extracted, not just those on the line render path, or whose center point is within the polygon. Not relevant for points}
  \item{method}{character. method for extracting values with points. The default is "simple", the alternative is "bilinear"}
  \item{...}{additional arguments. None implemented}  
}

\value{
matrix if \code{y} is a \code{SpatVector}, otherwise a vector.
}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
r[c(1:25, 31:100)] <- NA
r <- ifel(r > 28, r + 10, r)

# get all non-NA cell numbers
cells(r)

# get cell numbers that match values
cells(r, c(28,38))

m <- cbind(x=c(0,10,-30), y=c(40,-10,20))
cellFromXY(r, m)

v <- vect(m)
cells(r, v)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
#cv <- cells(r, v) 

#z  <- cells(r,ext(v))
#xy <- xyFromCell(r, z)
#plot(v)
#points(xy)

}


\keyword{methods}
\keyword{spatial}
