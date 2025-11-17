\name{boundaries}

\alias{boundaries}
\alias{boundaries,SpatRaster-method}
  
\title{Detect boundaries (edges)}

\description{
Detect boundaries (edges). Boundaries are cells that have more than one class in the 4 or 8 cells surrounding it, or, if \code{classes=FALSE}, cells with values and cells with \code{NA}. 
}

\usage{
\S4method{boundaries}{SpatRaster}(x, classes=FALSE, inner=TRUE, directions=8,
		falseval=0, ignoreNA=FALSE, filename="", ...)
}

\arguments{
\item{x}{SpatRaster}
\item{inner}{logical. If \code{TRUE}, "inner" boundaries are returned, else "outer" boundaries are returned}
\item{classes}{character. Logical. If \code{TRUE} all different values are (after rounding) distinguished, as well as \code{NA}. If \code{FALSE} (the default) only edges between \code{NA} and non-\code{NA} cells are considered}
\item{directions}{integer. Which cells are considered adjacent? Should be 8 (Queen's case) or 4 (Rook's case)}
\item{falseval}{numeric. The value to use for cells that are not a boundary and not \code{NA}}
\item{ignoreNA}{logical. If \code{TRUE} and \code{classes=TRUE} external boundaries (with \code{NA} cells) are ignored, only boundaries between classes are returned (and the value of argument \code{inner} is irrelevant)}
\item{filename}{character. Output filename}
\item{...}{options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster. Cell values are either 1 (a boundary) or 0 (not a boundary), or \code{NA}
}

\seealso{ \code{\link{focal}}, \code{\link{patches}} }

\examples{
r <- rast(nrows=18, ncols=36, xmin=0)
r[150:250] <- 1
r[251:450] <- 2
bi <- boundaries(r)
bo <- boundaries(r, inner=FALSE)
bc <- boundaries(r, classes=TRUE)
#plot(bc)
}

\keyword{methods}
\keyword{spatial}
