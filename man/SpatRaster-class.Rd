\name{SpatRaster-class}

\docType{class}

\alias{SpatRaster}
\alias{SpatRaster-class}
\alias{Rcpp_SpatRaster-class}
\alias{PackedSpatRaster-class}

\alias{SpatRasterCollection}
\alias{SpatRasterDataset}
\alias{SpatRasterCollection-class}
\alias{SpatRasterDataset-class}

\alias{RasterSource}
\alias{RasterSource-class}
\alias{Rcpp_RasterSource-class}

\alias{SpatCategories}
\alias{SpatCategories-class}
\alias{Rcpp_SpatCategories-class}

\alias{show,SpatRaster-method}


\title{SpatRaster class}

\description{
A \code{SpatRaster} represents a rectangular part of the world that is sub-divided into rectangular cells of equal area (in terms of the units of the coordinate reference system). For each cell can have multiple values ("layers"). An object of the \code{SpatRaster} class can point to one or more files on disk that hold the cell values, and/or it can hold these values in memory. These objects can be created with the \code{\link{rast}} method. 

A \code{SpatRasterDataset} is a collection of sub-datasets, where each is a \code{SpatRaster} with the same raster geometry (extent and resolution) and coordinate reference system. Sub-datasets are often used to represent different variables (e.g. temperature and precipitation), or a fourth dimension (e.g. height, depth or time) if the sub-datasets already have three dimensions (multiple layers). It can be created with the \code{\link{sds}} method.

A \code{SpatRasterCollection} is a collection of SpatRasters with no restriction in the extent or other geometric parameters. It can be created with the \code{\link{sprc}} method.

}


\examples{
rast()
}

\keyword{classes}
\keyword{spatial}
