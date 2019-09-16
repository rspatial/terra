\name{slope}

\alias{slope}
\alias{slope,SpatRaster-method}

\title{Compute slopes}

\description{
Compute slopes from elevation data. The elevation values should be in map units (typically meter) for projected (planar) raster data. They should be in meters when the coordinate reference system (CRS) is longitude/latitude. 
}

\usage{
\S4method{slope}{SpatRaster}(x, neighbors=8, unit="degrees", filename="", overwrite=FALSE, wopt=list(), ...)  
}

\arguments{
  \item{x}{Single layer SpatRaster with elevation values. Values should have the same unit as the map units, or in meters when the crs is longitude/latitude}
  \item{unit}{Character. 'degrees' or 'tangent'}
  \item{neighbors}{Integer. Indicating how many neighboring cells to use to compute slope for any cell. Either 8 (queen case) or 4 (rook case)}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}


\details{
When \code{neighbors=4} slope and aspect are computed according to Fleming and Hoffer (1979) and Ritter (1987). When \code{neigbors=8}, slope and aspect are computed according to Horn (1981). The Horn algorithm may be best for rough surfaces, and the Fleming and Hoffer algorithm may be better for smoother surfaces (Jones, 1997; Burrough and McDonnell, 1998).

If slope = 0, aspect is set to 0.5*pi radians (or 90 degrees if unit="degrees"). When computing slope or aspect, the crs of SpatRaster \code{x} must be known (may not be \code{NA}), to be able to safely differentiate between planar and longitude/latitude data.
}


\references{
Burrough, P., and R.A. McDonnell, 1998. Principles of Geographical Information Systems. Oxford University Press.

Fleming, M.D. and Hoffer, R.M., 1979. Machine processing of landsat MSS data and DMA topographic data for forest cover type mapping. LARS Technical Report 062879. Laboratory for Applications of Remote Sensing, Purdue University, West Lafayette, Indiana.

Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69:14-47

Jones, K.H., 1998. A comparison of algorithms used to compute hill slope as a property of the DEM. Computers & Geosciences 24: 315-323 

Ritter, P., 1987. A vector-based slope and aspect generation algorithm. Photogrammetric Engineering and Remote Sensing 53: 1109-1111
}


\examples{
\dontrun{
library(geodata)
elv <- elevation_30s("Switzerland", path=".")
x <- slope(elv, unit="degrees")
plot(x)
}
}

\keyword{spatial}
