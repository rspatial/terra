\name{writeCDF}

\alias{writeCDF}
\alias{writeCDF,SpatRasterDataset-method}
\alias{writeCDF,SpatRaster-method}

\title{Write raster data to a NetCDF file}

\description{
Write a SpatRaster or SpatRasterDataset to a NetCDF file. 

When using a SpatRasterDataset, the varname, longname, and unit should be set in the object (see examples).
}

\usage{
\S4method{writeCDF}{SpatRaster}(x, filename, varname, longname="", unit="", overwrite=FALSE,
 datatype="double", NAflag=-9999, ...)

\S4method{writeCDF}{SpatRasterDataset}(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{filename}{character. Output filename}
  \item{varname}{character}
  \item{longname}{character}
  \item{unit}{character}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{datatype}{character. One of "double", "float", "integer", "short", "byte" or "char"}
  \item{NAflag}{numeric}
  \item{...}{additional arguments passed on to \code{\link[ncdf4]{ncvar_def}}}
}

\value{
SpatRaster (invisibly). This method is used for the side-effect of writing values to a file.
}

\seealso{
see \code{\link{writeRaster}} for other file formats, or for writing netCDF with the GDAL library.
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
fname <- paste0(tempfile(), ".nc")

x <- writeCDF(r, fname, overwrite=TRUE, varname="alt", longname="elevation in m above sea level", unit="m")

a <- rast(ncol=5, nrow=5, nl=50)
values(a) <- 1:prod(dim(a))
time(a) <- as.Date("2020-12-31") + 1:nlyr(a)
x <- writeCDF(a, fname, overwrite=TRUE, varname="power", longname="my nice data", unit="U/Pa")

# for CRAN
file.remove(fname)
}

\keyword{ spatial }
\keyword{ methods }

