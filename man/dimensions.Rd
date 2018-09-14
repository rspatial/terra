\name{dimensions}

\docType{methods}

\alias{ncol}
\alias{nrow}
\alias{ncell}
\alias{nlayer}
\alias{dim}
\alias{ncol<-}
\alias{nrow<-}
\alias{ncell,ANY-method}
\alias{nrow,SpatRaster-method}
\alias{ncol,SpatRaster-method}
\alias{nlayer,SpatRaster-method}
\alias{ncell,SpatRaster-method}
\alias{length,SpatRaster-method}
\alias{nrow<-,SpatRaster-method}
\alias{ncol<-,SpatRaster-method}
\alias{dim,SpatRaster-method}
\alias{dim<-,SpatRaster-method}
\alias{length,SpatRaster-method}

\alias{res}
\alias{xres}
\alias{yres}
\alias{res,SpatRaster-method}
\alias{xres,SpatRaster-method}
\alias{yres,SpatRaster-method}


\title{Dimensions of a SpatRaster object}

\description{
Get or set the dimensions (number of rows, columns, cells, layers, length) of a SpatRaster object.
When setting dimensions, raster values are lost.
}

\usage{
\S4method{ncol}{SpatRaster}(x)
\S4method{nrow}{SpatRaster}(x)
\S4method{ncell}{SpatRaster}(x)
\S4method{nlayer}{SpatRaster}(x)
\S4method{res}{SpatRaster}(x)
\S4method{xres}{SpatRaster}(x)
\S4method{yres}{SpatRaster}(x)

\S4method{ncol}{SpatRaster}(x)<-value
\S4method{nrow}{SpatRaster}(x)<-value

}


\arguments{
\item{x}{SpatRaster object}
\item{value}{For ncol and nrow: postive integer. For res: one or two postive numbers }  
}

\value{
integer
}

\seealso{ \link{ext}}


\examples{
r <- georst()
ncell(r)
ncol(r)
nrow(r)
dim(r)

# nrow(r) <- 18
# ncol(r) <- 36
# equivalent to
# dim(r) <- c(18, 36) 

dim(r)
#dim(r) <- c(10, 10, 5)
dim(r)

xres(r)
yres(r)
res(r)

#res(r) <- 1/120
# set yres differently
#res(r) <- c(1/120, 1/60)
}

\keyword{spatial}

