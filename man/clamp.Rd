\name{clamp}

\alias{clamp}
\alias{clamp,SpatRaster-method}

\title{Clamp values}

\description{
Clamp values to a minimum and maximum value. That is, all values below a lower threshold value and above the upper threshold value become either \code{NA}, or, if \code{values=TRUE}, become the threshold value  
}

\usage{
\S4method{clamp}{SpatRaster}(x, lower=-Inf, upper=Inf, values=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{lower}{numeric. lowest value}
  \item{upper}{numeric. highest value}
  \item{values}{logical. If \code{FALSE} values outside the clamping range become \code{NA}, if \code{TRUE}, they get the extreme values}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\seealso{
\code{\link{classify}}
}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
rc <- clamp(r, 25, 75) 
rc
}

\keyword{spatial}

