\name{rangeFill}

\alias{rangeFill}
\alias{rangeFill,SpatRaster-method}

\title{Fill layers with a range}

\description{
Fill layers with cell-varying ranges defined by a start and end SpatRaster. The range must start at 1 and end at a user-defined maximum. Output values are either zero (not in the range) or one (in the range).

For example, for a cell with \code{start=3}, \code{end=5} and with \code{limit=8}, the output for that cell would be \code{0,0,1,1,1,0,0,0}
}

\usage{
\S4method{rangeFill}{SpatRaster}(x, limit, circular=FALSE, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster with at two layers. The cell values of the first layer indicate the start of the range (1 based); the cell values are indicate the end of the range}

\item{limit}{numeric > 1. The range size}

\item{circular}{logical. If \code{TRUE} the values are considered circular, such as the days of the year. In that case, if first > last the layers used are c(first:limit, 1:last). Otherwise, if \code{circular=FALSE}, such a range would be considered invalid and \code{NA} would be used} 

\item{filename}{character. Output filename}

\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}

}

\value{SpatRaster}


\seealso{ \code{\link{rapp}} } 


\examples{

x <- y <- rast(ncol=2, nrow=2)
values(x) <- c(NA, 1:3)
values(y) <- c(NA, 4:6)

r <- rangeFill(c(x, y), 8)
}


\keyword{spatial}

