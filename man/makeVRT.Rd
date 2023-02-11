\name{makeVRT}

\alias{makeVRT}

\title{Make a VRT header file}

\description{ 
Create a VRT header file for a "flat binary" raster file that needs a header file to be able to read it, but does not have it. 
}

\usage{
makeVRT(filename, nrow, ncol, nlyr=1, extent, xmin, ymin, xres, yres=xres, xycenter=TRUE,
   crs="+proj=longlat", lyrnms="", datatype, NAflag=NA, bandorder="BIL", byteorder="LSB",
   toptobottom=TRUE, offset=0, scale=1)
}

\arguments{
  \item{filename}{character. raster filename (without the ".vrt" extension)}
  \item{nrow}{positive integer, the number of rows}
  \item{ncol}{positive integer, the number of columns}
  \item{nlyr}{positive integer, the number of layers}
  \item{extent}{SpatExtent or missing}
  \item{xmin}{numeric. minimum x coordinate (only used if \code{extent} is missing)}
  \item{ymin}{numeric. minimum y coordinate (only used if \code{extent} is missing)}
  \item{xres}{positive number. x resolution}
  \item{yres}{positive number. y resolution)}
  \item{xycenter}{logical. If \code{TRUE}, \code{xmin} and \code{xmax} represent the coordinates of the center of the extreme cell, in stead of the coordinates of the outside corner. Only used of \code{extent} is missing}
  \item{crs}{character. Coordinate reference system description}
  \item{lyrnms}{character. Layer names}
  \item{datatype}{character. One of "INT2S", "INT4S", "INT1U", "INT2U", "INT4U", "FLT4S", "FLT8S". If missing, this is guessed from the file size (INT1U for 1 byte per value, INT2S for 2 bytes and FLT4S for 4 bytes per value). This may be wrong because, for example, 2 bytes per value may in fact be INT2U (with the U for unsigned) values}
  \item{NAflag}{numeric. The value used as the "NA flag"}
  \item{bandorder}{character. One of "BIL", "BIP", or "BSQ". That is Band Interleaved by Line, or by Pixel, or Band SeQuential}
  \item{byteorder}{character. One of "LSB", "MSB". "MSB" is common for files generated on Linux systems, whereas "LSB" is common for files generated on windows}
  \item{toptobottom}{logical. If \code{FALSE}, the values are read bottom to top}
  \item{offset}{numeric. offset to be applied}
  \item{scale}{numeric. scale to be applied}
}

\value{
character (.VRT filename)
}

\seealso{\code{\link{vrt}} to create a vrt for a collection of raster tiles}


\keyword{methods}
\keyword{spatial}
