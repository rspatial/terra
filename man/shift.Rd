\name{shift}

\docType{methods}

\alias{shift}
\alias{shift,SpatRaster-method}
\alias{shift,SpatExtent-method}
\alias{shift,SpatVector-method}
  
\title{Shift}

\description{
Shift the a SpatRaster, SpatVector or SpatExtent to another location.
}

\usage{
\S4method{shift}{SpatRaster}(x, dx=0, dy=0, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{shift}{SpatVector}(x, dx=0, dy=0, ...)

\S4method{shift}{SpatExtent}(x, dx=0, dy=0, ...)
}

\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{dx}{numeric. The shift in horizontal direction}
  \item{dy}{numeric. The shift in vertical direction}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}
 
\value{
Same as \code{x}
}


\seealso{ \code{\link{flip}}, \code{\link{rotate}}}

\examples{
r <- rast(xmn=0, ymn=0, xmx=1, ymx=1)
r <- shift(r, dx=1, dy=-1)

e <- ext(r)
shift(e, 5, 5)
}

\keyword{spatial}
