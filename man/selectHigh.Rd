\name{selectHighest}

\docType{methods}

\alias{selectHighest}
\alias{selectHighest,SpatRaster-method}

\title{select cells with high or low values}

\description{
Identify the cells that have the \code{n} highest or lowest values in the first layer of a SpatRaster.

With \code{ties="all"} (the default), all cells that have one of these values are selected, so the number of \code{TRUE} cells can be larger than \code{n} when there are ties (as is common with integer data). With \code{ties="first"}, exactly \code{n} cells are selected; ties are broken by choosing the smallest cell numbers. 
}

\usage{
\S4method{selectHighest}{SpatRaster}(x, n, low=FALSE, ties="first", filename="", ...)
}

\arguments{
  \item{x}{SpatRaster. Only the first layer is processed}
  \item{n}{positive integer. The number of highest or lowest values to select}
  \item{low}{logical. If \code{TRUE}, the lowest values are selected instead of the highest values}
  \item{ties}{character. Either \code{"all"} to select all cells that have one of the \code{n} highest or lowest values (the number of cells can be larger than \code{n} when there are ties), or \code{"first"} to select exactly \code{n} cells, breaking ties by choosing the smallest cell numbers}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}

\seealso{\code{\link{where.max}}, \code{\link{where.min}}, \code{\link{ifel}}, \code{\link{selectSize}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- selectHighest(r, 1000)
y <- selectHighest(r, 1000, TRUE)

m <- merge(y-1, x)
levels(m) <- data.frame(id=0:1, elevation=c("low", "high"))
plot(m)
}

\keyword{spatial}
