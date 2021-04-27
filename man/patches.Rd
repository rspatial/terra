\name{patches}

\alias{patches}
\alias{patches,SpatRaster-method}
  
\title{Detect patches (clumps) of cells}

\description{
Detect patches (clumps). Patches are groups of cells that are surrounded by cells that are \code{NA}. Set \code{zeroAsNA} to \code{TRUE} to also identify patches separated by cells with values of zero. 
}

\usage{
\S4method{patches}{SpatRaster}(x, directions=4, zeroAsNA=FALSE, filename="", ...)
}

\arguments{
\item{x}{SpatRaster}
\item{directions}{integer indicating which cells are considered adjacent. Should be 8 (Queen's case) or 4 (Rook's case)}
  \item{zeroAsNA}{logical. If \code{TRUE} treat cells that are zero as if they were \code{NA}}
  \item{filename}{character. Output filename}
  \item{...}{options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster. Cell values are either a patch number
}

\seealso{ \code{\link{focal}}, \code{\link{boundaries}} }

\examples{
r <- rast(nrow=18, ncol=36, xmin=0)
r[1:2, 5:8] <- 1
r[5:8, 2:6] <- 1
r[7:12, 22:36] <- 1
r[15:16, 18:29] <- 1
p <- patches(r)

r <- rast(nrow=10, ncol=10, xmin=0)
r[] <- 0
r[3, 3] <- 10
r[4, 4] <- 10
r[5, 5:8] <- 12
r[6, 6:9] <- 12

# remove zeros if need be with zeroAsNA

p4 <- patches(r, zeroAsNA=TRUE)
p8 <- patches(r, 8, zeroAsNA=TRUE)

# patches for different values
# remove zeros
rr <- classify(r, cbind(0, NA))
# first make layers for each value
s <- segregate(rr, keep=TRUE, other=NA)
p <- patches(s)

}

\keyword{methods}
\keyword{spatial}
