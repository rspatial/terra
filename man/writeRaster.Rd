\name{writeRaster}

\alias{writeRaster,SpatRaster,character-method}

\alias{writeRaster}

\title{Write raster data to a file}

\description{
Write a SpatRaster object to a file. 
}

\usage{
\S4method{writeRaster}{SpatRaster,character}(x, filename, overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename. Can be a single filename, or as many filenames as \code{nlyr(x)} to write a file for each layer}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments for for writing files. See Details }
}

\value{
SpatRaster. This function is used for the side-effect of writing values to a file.
}

\seealso{
see \code{\link{writeCDF}} for writing NetCDF files.
}

\details{
In writeRaster, and in other methods that generate SpatRaster objects, options for writing raster files to disk can be provided as additional arguments or, in a few cases, as the \code{wopt} argument (a named list) if the additional arguments are already used for a different purpose. See \code{\link{terraOptions}} to get or set default values. The following options are available: 

\tabular{ll}{
\bold{name} \tab \bold{description}\cr

\code{datatype}\tab values accepted are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". With GDAL >= 3.5 you can also use "INT8U" and "INT8S". And with GDAL >= 3.7 you can use also use "INT1S". See \code{\link{gdal}} to discover the GDAL version you are using.
The first three letters indicate whether the datatype is an integer (whole numbers) of a real number ("float", decimal numbers), the fourth character indicates the number of bytes used for each number. Higher values allow for storing larger numbers and/or more precision; but create larger files. The  "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (zero and positive values only).\cr

\code{filetype}\tab file format expresses as \href{https://gdal.org/drivers/raster/index.html}{GDAL driver names}. If this argument is not supplied, the driver is derived from the filename. You can use \code{gdal(drivers=TRUE)} to see what drivers are available in your installation\cr

\code{gdal}\tab GDAL driver specific datasource creation options. See the GDAL documentation. For example, with the \href{https://gdal.org/drivers/raster/gtiff.html}{GeoTiff file format} you can use \code{gdal=c("COMPRESS=DEFLATE", "TFW=YES")}.\cr

\code{tempdir}\tab the path where temporary files are to be written to.\cr

\code{progress}\tab positive integer. If the number of chunks is larger, a progress bar is shown.\cr

\code{memfrac}\tab numeric between 0 and 0.9 (higher values give a warning). The fraction of available RAM that terra is allowed to use.\cr

\code{memmax}\tab memmax - the maximum amount of RAM (in GB) that terra can use when processing a raster dataset. Should be less than what is detected (see \code{\link{mem_info}}, and higher values are ignored. Set it to a negative number or NA to ignore this value). \cr

\code{names}\tab output layer names.\cr

\code{NAflag}\tab numeric. value to represent missing (\code{NA} or \code{NaN}) values. See note\cr

\code{scale}\tab numeric. Cell values written to disk are divided by this value (default is 1). See \code{\link{scoff}}\cr

\code{offset}\tab numeric. Value that is subtracted from the cell values written to disk (default is 0). See 
\code{\link{scoff}} \cr

\code{verbose}\tab logical. If \code{TRUE} debugging information is printed\cr

\code{steps}\tab positive integers. In how many steps (chunks) do you want to process the data (for debugging)\cr

\code{todisk}\tab logical. If \code{TRUE} processing operates as if the dataset is very large and needs to be written to a temporary file (for debugging).\cr

}
}


\note{
GeoTiff files are, by default, written with LZW compression. If you do not want compression, use \code{gdal="COMPRESS=NONE"}.

When writing integer values the lowest available value (given the datatype) is used to represent \code{NA} for signed types, and the highest value is used for unsigned values. This can be a problem with byte data (between 0 and 255) as the value 255 is reserved for \code{NA}. To keep the value 255, you need to set another value as \code{NAflag}, or do not set a \code{NAflag} (with \code{NAflag=NA})
}

\examples{
r <- rast(nrows=5, ncols=5, vals=1:25)

# create a temporary filename for the example
f <- file.path(tempdir(), "test.tif")

writeRaster(r, f, overwrite=TRUE)

writeRaster(r, f, overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"), datatype='INT1U')

## Or with a wopt argument:

writeRaster(r, f, overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='INT1U'))

## remove the file
unlink(f)
}

\keyword{ spatial }
\keyword{ methods }

