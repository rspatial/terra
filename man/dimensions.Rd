\name{dimensions}

\docType{methods}

\alias{ncol}
\alias{nrow}
\alias{ncell}
\alias{nlyr}
\alias{nsrc}
\alias{dim}
\alias{ncol<-}
\alias{nrow<-}
\alias{ncell,ANY-method}
\alias{nrow,SpatRaster-method}
\alias{ncol,SpatRaster-method}
\alias{nrow<-,SpatRaster,numeric-method}
\alias{ncol<-,SpatRaster,numeric-method}

\alias{nlyr,SpatRaster-method}
\alias{nsrc,SpatRaster-method}
\alias{ncell,SpatRaster-method}
\alias{size,SpatRaster-method}
\alias{dim,SpatRaster-method}
\alias{dim<-,SpatRaster-method}
\alias{size,SpatVector-method}

\alias{res}
\alias{xres}
\alias{yres}
\alias{res,SpatRaster-method}
\alias{xres,SpatRaster-method}
\alias{yres,SpatRaster-method}

\alias{nrow,SpatVector-method}
\alias{ncol,SpatVector-method}

\title{Dimensions of a SpatRaster or SpatVector}

\description{
Get the number of rows (nrow), columns (ncol), cells (ncell), layers (nlyr), sources (nsrc) or the spatial resolution of a SpatRaster.
Set the number of rows or columns or layers. When setting dimensions, all sources and values are dropped.
}

\usage{
\S4method{ncol}{SpatRaster}(x)
\S4method{nrow}{SpatRaster}(x)
\S4method{nlyr}{SpatRaster}(x)
\S4method{ncell}{SpatRaster}(x)
\S4method{size}{SpatRaster}(x)
\S4method{nsrc}{SpatRaster}(x)
\S4method{res}{SpatRaster}(x)

ncol(x, ...) <- value
nrow(x, ...) <- value

\S4method{xres}{SpatRaster}(x)
\S4method{yres}{SpatRaster}(x)

\S4method{ncol}{SpatVector}(x)
\S4method{nrow}{SpatVector}(x)
\S4method{size}{SpatVector}(x, ...)

}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{value}{For ncol and nrow: positive integer. For res: one or two positive numbers }  
\item{...}{additional arguments. None implemented}
}

\value{
integer
}

\seealso{ \link{ext}}


\examples{
r <- rast()
ncell(r)
ncol(r)
nrow(r)
dim(r)
nsrc(r)

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

