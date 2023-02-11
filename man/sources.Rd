\name{sources}

\docType{methods}

\alias{sources}
\alias{sources,SpatRaster-method}
\alias{sources,SpatRasterCollection-method}
\alias{sources,SpatRasterDataset-method}
\alias{sources,SpatVector-method}
\alias{sources,SpatVectorProxy-method}

\alias{hasValues}
\alias{hasValues,SpatRaster-method}
\alias{inMemory}
\alias{inMemory,SpatRaster-method}


\title{Data sources of a SpatRaster}

\description{
Get the data sources of a SpatRaster or SpatVector or related object. Sources are either files (or similar resources) or "", meaning that they are in memory. You can use \code{hasValues} to check if in-memory layers actually have cell values. 
}

\usage{
\S4method{sources}{SpatRaster}(x, nlyr=FALSE, bands=FALSE)

\S4method{sources}{SpatVector}(x)

\S4method{hasValues}{SpatRaster}(x)

\S4method{inMemory}{SpatRaster}(x, bylayer=FALSE)
}

\arguments{
\item{x}{SpatRaster, SpatRasterCollection, SpatVector or SpatVectorProxy}
\item{nlyr}{logical. If \code{TRUE} for each source, the number of layers is returned}
\item{bands}{logical. If \code{TRUE} for each source, the "bands" used, that is, the layer number in the source file, are returned}
\item{bylayer}{logical. If \code{TRUE} a value is returned for each layer instead of for each source}
}

\value{
A vector of filenames, or \code{""} when there is no filename, if \code{nlyr} and \code{bands} are both \code{FALSE}. Otherwise a \code{data.frame}
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
s <- rast(r)
values(s) <- 1:ncell(s)
rs <- c(r,r,s,r)
sources(rs)
hasValues(r)
x <- rast()
hasValues(x)
}


\keyword{methods}
\keyword{spatial}

