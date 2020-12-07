\name{writeCDF}

\alias{writeCDF}
\alias{writeCDF,SpatDataSet-method}
\alias{writeCDF,SpatRaster-method}


\title{Write raster data to a NetCDF file}

\description{
Write a SpatRaster or SpatDataSet to a NetCDF file.
}

\usage{
\S4method{writeCDF}{SpatRaster}(x, filename,  varname, longname="", units="", overwrite=FALSE, datatype="double", NAflag=-9999, ...)

\S4method{writeCDF}{SpatDataSet}(x, filename, overwrite=FALSE, datatype="double", NAflag=-9999, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatDataSet}
  \item{filename}{character. Output filename}
  \item{varname}{character)
  \item{longname}{character)
  \item{units}{character)
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{datatype}{character. One of "double", "float", "integer", "short", "byte" or "char"}
  \item{NAflag}{numeric}
  \item{...}{additional arguments passed on to \code{\link[ncdf4]{ncvar_def}}}
}

\value{
SpatRaster (invisibly). This method is used for the side-effect of writing values to a file.
}

\seealso{
see \code{\link{writeRaster}} for other file formats.
}


\keyword{ spatial }
\keyword{ methods }

