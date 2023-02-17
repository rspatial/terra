\name{rapp}

\docType{methods}

\alias{rapp}
\alias{rapp,SpatRaster-method}

\title{Range-apply}

\description{
Apply a function to a range of the layers of a SpatRaster that varies by cell. The range is specified for each cell with one or two SpatRasters (arguments \code{first} and \code{last}). For either \code{first} or \code{last} you can use a single number instead. 

You cannot use single numbers for both \code{first} and \code{last} because in that case you could use \code{\link{app}} or \code{\link{Summary-methods}}, perhaps \code{\link{subset}}ting the layers of a SpatRaster.

See \code{\link{selectRange}} to create a new SpatRaster by extracting one or more values starting at a cell-varying layer.
}

\usage{
\S4method{rapp}{SpatRaster}(x, first, last, fun, ..., allyrs=FALSE, fill=NA, 
        clamp=FALSE, circular=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{first}{SpatRaster or positive integer between 1 and nlyr(x), indicating the first layer in the range of layers to be considered}
  \item{last}{SpatRaster or positive integer between 1 and nlyr(x), indicating the last layer in the range to be considered}  
  \item{fun}{function to be applied}
  \item{...}{additional arguments passed to \code{fun}} 
  \item{allyrs}{logical. If \code{TRUE}, values for all layers are passed to \code{fun} but the values outside of the range are set to \code{fill}}
  \item{fill}{numeric. The fill value for the values outside of the range, for when \code{allyrs=TRUE}} 
  \item{clamp}{logical. If \code{FALSE} and the specified range is outside \code{1:nlyr(x)} all cells are considered \code{NA}. Otherwise, the invalid part of the range is ignored} 
  \item{circular}{logical. If \code{TRUE} the values are considered circular, such as the days of the year. In that case, if first > last the layers used are c(first:nlyr(x), 1:last). Otherwise, the range would be considered invalid and \code{NA} would be returned} 
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{selectRange}}, \code{\link{app}}, \code{\link{Summary-methods}}, \code{\link{lapp}}, \code{\link{tapp}}}

\examples{
r <- rast(ncols=9, nrows=9)
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

