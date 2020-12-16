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
\S4method{writeCDF}{SpatRaster}(x, filename, varname, longname="", unit="", ...)

\S4method{writeCDF}{SpatRasterDataset}(x, filename, overwrite=FALSE, zname="time", missval=-9999, prec="float", compression=NA, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{filename}{character. Output filename}
  \item{varname}{character. Name of the dataset}
  \item{longname}{character. Long name of the dataset}
  \item{unit}{character. Unit of the data}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{zname}{character. The name of the "time" dimension}
  \item{missval}{numeric, the number used to indicate missing values}
  \item{prec}{character. One of "double", "float", "integer", "short", "byte" or "char"}
  \item{compression}{see \code{\link[ncdf4]{ncvar_def}}}
  \item{...}{additional arguments passed on to \code{\link[ncdf4]{ncvar_def}}}
}

\value{
logical (invisibly) indicating success
}

\seealso{
see \code{\link{writeRaster}} for writing other file formats
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
fname <- paste0(tempfile(), ".nc")

writeCDF(r, fname, overwrite=TRUE, varname="alt", longname="elevation in m above sea level", unit="m")

a <- rast(ncol=5, nrow=5, nl=50)
values(a) <- 1:prod(dim(a))
time(a) <- as.Date("2020-12-31") + 1:nlyr(a)
writeCDF(a, fname, overwrite=TRUE, varname="power", longname="my nice data", unit="U/Pa")


b <- sqrt(a)
s <- sds(a, b)
names(s) <- c("temp", "prec")
longnames(s) <- c("temperature (C)", "precicpiation (mm)")
units(s) <- c("C", "mm")
writeCDF(s, fname, overwrite=TRUE)

# for CRAN
file.remove(fname)
}

\keyword{ spatial }
\keyword{ methods }

