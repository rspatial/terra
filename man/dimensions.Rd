\name{dimensions}

\docType{methods}

\alias{size}
\alias{length}
\alias{ncol}
\alias{nrow}
\alias{ncell}
\alias{nlyr}
\alias{nsrc}
\alias{dim}
\alias{res}
\alias{xres}
\alias{yres}
\alias{ncol<-}
\alias{nrow<-}
\alias{nlyr<-}
\alias{res<-}
\alias{size,SpatRaster-method}
\alias{ncell,ANY-method}
\alias{nrow,SpatRaster-method}
\alias{ncol,SpatRaster-method}
\alias{nrow,SpatRasterDataset-method}
\alias{ncol,SpatRasterDataset-method}
\alias{nrow,SpatRasterCollection-method}
\alias{ncol,SpatRasterCollection-method}
\alias{nrow<-,SpatRaster,numeric-method}
\alias{ncol<-,SpatRaster,numeric-method}
\alias{nlyr,SpatRasterDataset-method}
\alias{res<-,SpatRaster,numeric-method}

\alias{nlyr,SpatRaster-method}
\alias{nlyr<-,SpatRaster,numeric-method}
\alias{nsrc,SpatRaster-method}
\alias{ncell,SpatRaster-method}
\alias{ncell,SpatRasterDataset-method}
\alias{length,SpatRasterDataset-method}
\alias{length,SpatRasterCollection-method}
\alias{nlyr,SpatRasterCollection-method}
\alias{dim,SpatRaster-method}
\alias{dim,SpatRasterDataset-method}
\alias{dim,SpatRasterCollection-method}
\alias{dim<-,SpatRaster-method}
\alias{dim,SpatVector-method}
\alias{dim,SpatVectorProxy-method}

\alias{res,SpatRasterDataset-method}
\alias{res,SpatRaster-method}
\alias{res<-,SpatRaster-method}
\alias{xres,SpatRaster-method}
\alias{yres,SpatRaster-method}

\alias{nrow,SpatVector-method}
\alias{ncol,SpatVector-method}

\alias{length,SpatVector-method}
\alias{length,SpatVectorCollection-method}


\title{Dimensions of a SpatRaster or SpatVector and related objects}

\description{
Get the number of rows (\code{nrow}), columns (\code{ncol}), cells (\code{ncell}), layers (\code{nlyr}), sources (\code{nsrc}), the size \code{size} (\code{nlyr(x)*ncell(x)}), or spatial resolution of a SpatRaster.

\code{length} returns the number of sub-datasets in a SpatRasterDataset or SpatVectorCollection.

For a SpatVector \code{length(x)} is the same as \code{nrow(x)}.

You can also set the number of rows or columns or layers. When setting dimensions, all cell values are dropped.
}

\usage{
\S4method{ncol}{SpatRaster}(x)

\S4method{nrow}{SpatRaster}(x)

\S4method{nlyr}{SpatRaster}(x)

\S4method{ncell}{SpatRaster}(x)

\S4method{nsrc}{SpatRaster}(x)

\S4method{ncol}{SpatRaster,numeric}(x)<-value

\S4method{nrow}{SpatRaster,numeric}(x)<-value

\S4method{nlyr}{SpatRaster,numeric}(x)<-value

\S4method{res}{SpatRaster}(x)

\S4method{res}{SpatRaster,numeric}(x)<-value

\S4method{xres}{SpatRaster}(x)

\S4method{yres}{SpatRaster}(x)

\S4method{ncol}{SpatVector}(x)

\S4method{nrow}{SpatVector}(x)

\S4method{length}{SpatVector}(x)

}


\arguments{
\item{x}{SpatRaster or SpatVector or related objects}
\item{value}{For ncol and nrow: positive integer. For res: one or two positive numbers }  
}

\value{
integer
}

\seealso{ \link{ext}}

\examples{
r <- rast()
ncol(r)
nrow(r)
nlyr(r)
dim(r)
nsrc(r)
ncell(r)

rr  <- c(r,r)
nlyr(rr)
nsrc(rr)
ncell(rr)

nrow(r) <- 18
ncol(r) <- 36
# equivalent to
dim(r) <- c(18, 36) 

dim(r)
dim(r) <- c(10, 10, 5)
dim(r)

xres(r)
yres(r)
res(r)

res(r) <- 1/120
# different xres and yres
res(r) <- c(1/120, 1/60)
}

\keyword{spatial}

