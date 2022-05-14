\name{project}

\alias{project}
\alias{project,SpatVector-method}
\alias{project,SpatRaster-method}
\alias{project,matrix-method}

\title{Change the coordinate reference system}

\description{
Change the coordinate reference system ("project") of a SpatVector, SpatRaster or a matrix with coordinates. 
}

\usage{
\S4method{project}{SpatVector}(x, y)

\S4method{project}{SpatRaster}(x, y, method, mask=FALSE, align=FALSE, 
		gdal=TRUE, res=NULL, origin=NULL, filename="", ...)

\S4method{project}{matrix}(x, from, to)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{if (\code{x} is a SpatRaster, the preferred approach is for \code{y} to be a SpatRaster as well, serving as a template for the geometry (extent and resolution) of the output SpatRaster. Alternatively, you can provide a coordinate reference system (CRS) description. 
    
  You can use the following formats to define coordinate reference systems: WKT,  PROJ.4 (e.g., \code{+proj=longlat +datum=WGS84}), or an EPSG code (e.g., \code{"epsg:4326"}). But note that the PROJ.4 notation has been deprecated, and you can only use if with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored. 
  
  If \code{x} is a SpatVector, you can provide a crs definition as discussed above, or any other object from which such a crs can be extracted with \code{\link{crs}}}
  
  \item{method}{character. Method used for estimating the new cell values of a SpatRaster. One of: 
  
	\code{near}: nearest neighbor. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for continuous values. This is used by default if the first layer of \code{x} is categorical.

	\code{bilinear}: bilinear interpolation. This is the default if the first layer of \code{x} is numeric (not categorical).
	
	\code{cubic}: cubic interpolation.
	
	\code{cubicspline}: cubic spline interpolation.
  }

  \item{mask}{logical. If \code{TRUE}, mask out areas outside the input extent (see example with Robinson projection)}
  \item{align}{logical. If \code{TRUE}, and \code{y} is a SpatRaster, the template is used for the spatial resolution and origin, but the extent is set such that all of the extent of \code{x} is included}

  \item{gdal}{logical. If \code{TRUE} the GDAL-warp algorithm is used. Otherwise a slower internal algorithm is used that may be more accurate if there is much variation in the cell sizes of the output raster. Only the \code{near} and \code{bilinear} algorithms are available for the internal algorithm}

  \item{res}{numeric. Can be used to set the resolution of the output raster if \code{y} is a CRS}

  \item{origin}{numeric. Can be used to set the origon of the output raster if \code{y} is a CRS}
  
  \item{filename}{character. Output filename}
  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}

  \item{from}{character. Coordinate reference system for \code{x}}

  \item{to}{character. Output coordinate reference system}
}


\value{
SpatVector or SpatRaster
}

\seealso{\code{\link{crs}}, \code{\link{resample}}}

\note{
User beware. Sadly, the PROJ.4 notation has been partly deprecated in the GDAL/PROJ library that is used by this function. You can still use this notation, but *only* with the the WGS84 datum. Other datums are silently ignored. 

When printing a Spat* object the PROJ.4 notation is shown because it is the most concise and clear format available. However, internally a WKT representation is used (see \code{\link{crs}}).

Transforming (projecting) raster data is fundamentally different from transforming vector data. Vector data can be transformed and back-transformed without loss in precision and without changes in the values. This is not the case with raster data. In each transformation the values for the new cells are estimated in some fashion. Therefore, if you need to match raster and vector data for analysis, you should generally transform the vector data. 

When using this method with a \code{SpatRaster}, the preferable approach is to provide a template \code{SpatRaster} as argument \code{y}. The template is then another raster dataset that you want your data to align with. If you do not have a template to begin with, you can do \code{project(x, crs)} and then manipulate the output to get the template you want. For example, where possible use whole numbers for the extent and resolution so that you do not have to worry about small differences in the future. You can use commands like \code{dim(z) = c(180, 360)} or \code{res(z) <- 100000}. 

The output resolution should be similar to the input resolution, but there is not "correct" resolution in raster transformation; but it is not obvious what this resolution is if you are using lon/lat data that spans a large North-South extent.

See \code{sf_proj_pipelines} to use specific "projection pipe lines"
}


\examples{
## SpatRaster
a <- rast(ncols=40, nrows=40, xmin=-110, xmax=-90, ymin=40, ymax=60, 
          crs="+proj=longlat +datum=WGS84")
values(a) <- 1:ncell(a)
newcrs="+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
b <- rast(ncols=94, nrows=124, xmin=-944881, xmax=935118, ymin=4664377, ymax=7144377, crs=newcrs)
w <- project(a, b)


## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
crs(v, proj=TRUE)
cat(crs(v), "\n")

crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"
p <- project(v, crs)
p

project(v, "EPSG:2169")
}

\keyword{spatial}

