\name{ar_info}

\alias{ar_info}

\title{Describe a multi-dimensional array file}

\description{
Describe a multi-dimensional array file (e.g., NetCDF)
}

\usage{
ar_info(x, what="describe", simplify=TRUE, filter=TRUE, array="")
}

\arguments{
  \item{x}{character. The name of a netcdf (or similar) raster file}
  \item{what}{character that (partially) matches "describe", "arrays" or "dimensions"}
  \item{simplify}{logical. If \code{TRUE} and \code{what="describe"}, simplify the output for readability}   
  \item{filter}{logical. If \code{TRUE} and \code{what="describe"} filter arrays that (probably) dimensions}   
  \item{array}{character. Required when \code{what="dimensions"}}   
}

\value{
character or data.frame
}

\seealso{\code{\link{describe}}}

\keyword{spatial}
