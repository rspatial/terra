\name{density}

\alias{density}
\alias{density,SpatRaster-method}

\docType{methods}

\title{Density plot}

\description{
Create density plots of the cell values of a SpatRaster
}

\usage{
\S4method{density}{SpatRaster}(x, maxcells=100000, plot=TRUE, main, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{maxcells}{the maximum number of (randomly sampled) cells to be used for creating the plot}
\item{plot}{if \code{TRUE} produce a plot, else return a density object}
\item{main}{character. Caption of plot(s)}
\item{...}{additional arguments passed to \code{\link{plot}}}
}

\value{
density plot (and a density object, returned invisibly if \code{plot=TRUE)} 
}

\examples{
logo <- rast(system.file("ex/logo.tif", package="terra"))
density(logo)
}

\keyword{spatial}
