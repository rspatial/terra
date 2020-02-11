\name{boundaries}

\alias{boundaries}
\alias{boundaries,SpatRaster-method}
  
\title{Detect boundaries (edges)}

\description{
Detect boundaries (edges). boundaries are cells that have more than one class in the 4 or 8 cells surrounding it, or, if \code{classes=FALSE}, cells with values and cells with \code{NA}. 
}

\usage{
\S4method{boundaries}{SpatRaster}(x, classes=FALSE, type="inner", 
         directions=8, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
\item{x}{SpatRaster}
\item{type}{character. "inner" or "outer"}
\item{classes}{character. Logical. If \code{TRUE} all different values are (after rounding) distinguished, as well as \code{NA}. If \code{FALSE} (the default) only edges between \code{NA} and non-\code{NA} cells are considered}
\item{directions}{integer. Which cells are considered adjacent? Should be 8 (Queen's case) or 4 (Rook's case)}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster. Cell values are either 1 (a border) or 0 (not a border), or \code{NA}
}

\seealso{ \code{\link{focal}}, \code{\link{clump}} }

\examples{
r <- rast(nrow=18, ncol=36, xmn=0)
v <- rep(NA, ncell(r))
v[150:250] <- 1
v[251:450] <- 2
values(r) <- v
plot( boundaries(r, type='inner') )
plot( boundaries(r, type='outer') )
plot( boundaries(r, classes=TRUE) )
}

\keyword{methods}
\keyword{spatial}
