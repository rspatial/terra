\name{project}

\alias{project}
\alias{project,SpatVector-method}
\alias{project,SpatRaster-method}

\title{Change the coordinate reference system}

\description{
Change the coordinate reference system ("project") of a SpatVector or SpatRaster. 
}

\usage{
\S4method{project}{SpatVector}(x, y)

\S4method{project}{SpatRaster}(x, y, method="bilinear", mask=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{if (\code{x} is a SpatRaster, the prefered approach is for \code{y} to be a SpatRaster as well, serving as a template for the output SpatRaster. Alternatively, you can provide a character variable containing a textual coordinate reference system (crs). 
    
  You can use the following formats to provide textual coordinate reference system (crs) definitions: WKT,  PROJ.4 (e.g., \code{+proj=longlat +datum=WGS84}), or an EPSG code (e.g., \code{"epsg:4326"}) to define the output coordinate reference system (crs). But note that when using PROJ.4 you now must use the WGS84 datum. Other datums are silently ignored. 
  
  If \code{x} is a SpatVector, you can provide a crs definition as discussed above, or any other object from which such a crs can be extracted with \code{\link{crs}}}
  
  \item{method}{character. Method used for estimating the new cell values. One of: 
  
	near: nearest neighbor. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for continuous values.
	bilinear: bilinear interpolation. Default.
	cubic: cubic interpolation.
	cubicspline: cubic spline interpolation.
  }

  \item{mask}{logical. If \code{TRUE}, mask out areas outside the input extent (see example with Robinson projection)}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatVector or SpatRaster
}

\seealso{\code{\link{crs}}, \code{\link{resample}}}

\note{
User beware. Sadly, the PROJ.4 notation has been partly deprecated in the GDAL/PROJ library that is used by this function. You can still use this notation, but *only* with the the WGS84 datum. Other datums are silently ignored. 

When printing a Spat* object the PROJ.4 notation is shown because it is the most concise and clear format available. However, internally a WKT representation is used (see \code{\link{crs}}).

Transforming (projecting) raster data is fundamentally different from transforming vector data. Vector data can be transformed and back-transformed without loss in precision and without changes in the values. This is not the case with raster data. In each transformation the values for the new cells are estiamted in some fashion. Therefore, if you need to match raster and vector data for analysis, you should generally transform the vector data. 

When using this method with a \code{SpatRaster}, the preferable approach is to provide a template \code{SpatRaster} as argument \code{y}. The template is then another raster dataset that you want your data to align with. If you do not have a template to begin with, you can do \code{project(x, crs)} and then manipulate the output to get the template you want. For example, where possible use whole numbers for the extent and resolution so that you do not have to worry about small differences in the future. You can use commands like \code{dim(z) = c(180, 360)} or \code{res(z) <- 100000}. 

The output resolution should be similar to the input resolution, but there is not "correct" resolution in raster transformation; but it is not obvious what this resolution is if you are using lon/lat data that spans a large North-South extent.
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
crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"
p <- project(v, crs)
p
}

\keyword{spatial}

