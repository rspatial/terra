\name{writeRaster}

\alias{writeRaster,SpatRaster,character-method}

\alias{writeRaster}

\title{Write raster data to a file}

\description{
Write a SpatRaster to a file. 
}

\usage{
\S4method{writeRaster}{SpatRaster,character}(x, filename, overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files. See Details }
  \item{...}{Additional arguments. None implemented}
}

\value{
SpatRaster (invisibly). This function is used for the side-effect of writing values to a file.
}

\details{
In many methods, options for writing raster files to disk can be provided with the \code{wopt} argument. \code{wopt} should be a named list. The following names are valid: \code{filetype} and \code{datatype}. 

Values for \code{filetype} are GDAL driver names. If this argument is not supplied, the driver is derived from the filename.

Values for \code{datatype} are "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". The first three letters indicate whether it is an integer (whole numbers) of float (decimal numbers), the fourth character indicates the number of bytes used (allowing for large numbers, and/or more precision), and the  "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (postive values only).
}

\keyword{ spatial }
\keyword{ methods }

