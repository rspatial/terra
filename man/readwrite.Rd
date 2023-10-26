\name{readwrite}

\docType{methods}

\alias{readStart}
\alias{readStop}
\alias{readValues}
\alias{writeStart}
\alias{writeStop}
\alias{writeValues}
\alias{fileBlocksize}
\alias{blocks}
\alias{readStart,SpatRaster-method}
\alias{readStop,SpatRaster-method}
\alias{readStart,SpatRasterDataset-method}
\alias{readStop,SpatRasterDataset-method}
\alias{readValues,SpatRaster-method}
\alias{writeStart,SpatRaster,character-method}
\alias{writeStop,SpatRaster-method}
\alias{writeValues,SpatRaster,vector-method}
\alias{blocks,SpatRaster-method}

\title{Read from, or write to, file}

\description{
Methods to read from or write chunks of values to or from a file. These are low level methods for programmers. Use writeRaster if you want to save an entire SpatRaster to file in one step. It is much easier to use. 

To write chunks, begin by opening a file with \code{writeStart}, then write values to it in chunks using the list that is returned by \code{writeStart}. When writing is done, close the file with \code{writeStop}. 

\code{blocks} only returns chunk size information. This can be useful when reading, but not writing, raster data. 
}

\usage{
\S4method{readStart}{SpatRaster}(x)

\S4method{readStop}{SpatRaster}(x)

\S4method{readValues}{SpatRaster}(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE, ...)

\S4method{writeStart}{SpatRaster,character}(x, filename="", overwrite=FALSE, n=4, sources="", ...)

\S4method{writeStop}{SpatRaster}(x)

\S4method{writeValues}{SpatRaster,vector}(x, v, start, nrows)

\S4method{blocks}{SpatRaster}(x, n=4)

fileBlocksize(x)
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename}

  \item{v}{vector with cell values to be written}
  \item{start}{integer. Row number (counting starts at 1) from where to start writing \code{v}}

 \item{row}{positive integer. Row number to start from, should be between 1 and nrow(x)}
 \item{nrows}{positive integer. How many rows?}
 \item{col}{positive integer. Column number to start from, should be between 1 and ncol(x)}
 \item{ncols}{positive integer. How many columns? Default is the number of columns left after the start column}

 \item{mat}{logical. If \code{TRUE}, values are returned as a numeric matrix instead of as a vector, except when \code{dataframe=TRUE}. If any of the layers of \code{x} is a factor, the level index is returned, not the label. Use \code{dataframe=TRUE} to get the labels}
 \item{dataframe}{logical. If \code{TRUE}, values are returned as a \code{data.frame} instead of as a vector (also if matrix is \code{TRUE})}
  
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{n}{positive integer indicating how many copies the data may be in memory at any point in time. This is used to determine how many blocks (large) datasets need to be read}
  \item{sources}{character. Filenames that may not be overwritten because they are used as input to the function. Can be obtained with \code{sources(x)}}

  \item{...}{
	For \code{writeStart}: additional arguments for writing files as in \code{\link{writeRaster}}

    For \code{readValues}: additional arguments for \code{\link{data.frame}} (and thus only relevant when \code{dataframe=TRUE})
  }
}

\value{
\code{readValues} returns a vector, matrix, or data.frame 

\code{writeStart} returns a list that can be used for processing the file in chunks.

The other methods invisibly return a logical value indicating whether they were successful or not. Their purpose is the side-effect of opening or closing files. 
}

\keyword{ spatial }
\keyword{ methods }
