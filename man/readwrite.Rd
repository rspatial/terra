\name{read and write}

\docType{methods}

\alias{readStart}
\alias{readStop}
\alias{writeStart}
\alias{writeStop}
\alias{writeValues}

\alias{readStart,SpatRaster-method}
\alias{readStop,SpatRaster-method}
\alias{writeStart,SpatRaster,character-method}
\alias{writeStop,SpatRaster-method}
\alias{writeStop,SpatRaster-method}
\alias{writeValues,SpatRaster,vector-method}

\title{read from or write to file}

\description{
Methods to read from or write chunks of values to or from a file. These are low level methods for programmers. Use writeRaster if you want to save an entire SpatRaster to file. 

Functions for writing chunks has to start at the first cell of a row (identified with argument \code{row}) and the values written must represent 1 or more entire rows. Begin by opening a file with \code{writeStart}, then write values to it in chunks. When writing is done close the file with \code{writeStop}. 

If you want to write all values of a Raster* object at once, you can also use \code{\link{writeRaster}} which is easier to use but more limited. 
}


\usage{
\S4method{readStart}{SpatRaster}(x, ...)
\S4method{readStop}{SpatRaster}(x)
\S4method{writeStart}{SpatRaster,character}(x, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...)
\S4method{writeStop}{SpatRaster}(x)
\S4method{writeValues}{SpatRaster,vector}(x, v, start)
}

\arguments{
  \item{x}{SpatRaster object}
  \item{filename}{Character. Output filename. Optional}
  \item{format}{Character. Output file format. See \code{\link{writeRaster}}}
  \item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{Additional arguments. None implemented}
  \item{v}{vector with cell values to be written}
  \item{start}{integer. Row number (counting starts at 1) from where to start writing \code{v}}

}


\value{
\code{writeStart} returns a list that can be used for processing the file in chunks.

The other methods invisibly return a logical value indicating whether they were succesful or not. Their purpose is the side-effect of opening or closing files. 
}

\keyword{ spatial }
\keyword{ methods }
