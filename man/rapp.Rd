\name{rapp}

\docType{methods}

\alias{rapp}
\alias{rapp,SpatRaster-method}

\title{Range-Apply}

\description{
Apply a function to a cell-varying range of the layers of a SpatRaster. The range is specified for each cell separately by two other SpatRasters (arguments \code{first} and \code{last}. For either \code{first} or \code{last} you can use a constant instead. 

If the specified range is outside \code{1:nlyr(x)}, \code{NA} is returned. You could use \code{\link{clamp}} to assure that this is the case.

See \code{\link{app}} or \code{\link{Summary-methods}} if you want to apply a function to all layers (perhaps after making a \code{\link{subset}}) of a SpatRaster.
}

\usage{
\S4method{rapp}{SpatRaster}(x, first, last, fun, ..., allyrs=FALSE, fill=NA, 
          filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{first}{SpatRaster or positive integer between 1 and nlyr(x), indicating the first layer in the range of layers to be considered}
  \item{last}{SpatRaster or positive integer between 1 and nlyr(x), indicating the last layer in the range to be considered}  
  \item{fun}{function to be applied}
  \item{...}{additional arguments passed to \code{fun}} 
  \item{allyrs}{logical. Should the values of all layers be passed to \code{fun}. The values outside of the range are set to \code{fill}}
  \item{fill}{numeric. The fill value for the the values outside of the range, for when \code{allyrs=TRUE}} 
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{app}}, \code{\link{Summary-methods}}, \code{\link{lapp}}, \code{\link{tapp}}}

\examples{
r <- rast(ncol=9, nrow=9)
values(r) <- 1:ncell(r)
s <- c(r, r, r, r, r, r)
s <- s * 1:6
s[1:2] <- NA
start <- end <- rast(r)
start[] <- 1:3
end[]   <- 4:6
a <- rapp(s, start, end, fun="mean")
b <- rapp(s, start, 2, fun="mean")

# cumsum from start to nlyr(x). return all layers
r <- rapp(s, start, nlyr(s), cumsum, allyrs=TRUE, fill=0)
# return only the final value
rr <- rapp(s, start, nlyr(s), function(i) max(cumsum(i)))
}

\keyword{methods}
\keyword{spatial}

