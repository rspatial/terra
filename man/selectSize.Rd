\name{selectSize}
\alias{selectSize}
\alias{selectSize,SpatVector-method}

\title{Select features by size}

\description{
Select the largest or smallest polygons or lines in a \code{SpatVector}, optionally providing a maximum number of selected features, filtering by a size threshold, dissolving adjacent features, and/or splitting multipart features into single parts.
}

\usage{
\S4method{selectSize}{SpatVector}(x, n=Inf, largest=TRUE, threshold=NULL,
           dissolve=FALSE, singleparts=TRUE)
}

\arguments{
  \item{x}{\code{SpatVector} of polygons or lines}
  \item{n}{optional integer value (default Inf) indicating the maximum number of largest or smallest features to select (if they pass the optional threshold)  }
  \item{largest}{if \code{TRUE} (default) select the largest feature(s), otherwise the smallest}
  \item{threshold}{optional size threshold (area in square meters for polygons, length in meters for lines) above which (if \code{largest=TRUE} or below which (if \code{largest=FALSE} features should be selected}
  \item{dissolve}{if \code{TRUE}, aggregate features before selection. Default \code{FALSE}}
  \item{singleparts}{if \code{TRUE} (default), split multipart features into individual parts before selection}
}

\details{
Polygon size is calculated using \code{\link{expanse}}; line size is calculated using \code{\link{perim}}. The threshold is applied before selecting the requested number of features. If fewer than \code{n} features satisfy the threshold, only qualifying features are returned.

If \code{dissolve=TRUE}, features are aggregated before size calculation. If \code{singleparts=TRUE}, multipart features are disaggregated before size calculation.
}

\value{\code{SpatVector} containing the selected features.}

\author{A. Márcia Barbosa}

\seealso{\code{\link{expanse}}, \code{\link{perim}}, \code{\link{selectHighest}}}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

# Select the five largest polygons
large5 <- selectSize(v, n=5)

# Select the three smallest polygons
small3 <- selectSize(v, n=3, largest=FALSE)

# Select up to 4 largest polygons satisfying a given minimum area
large <- selectSize(v, n=4, threshold=260e+6)

\dontrun{
plot(v)
polys(large5, col="red")
polys(small3, col="blue")
polys(large, col="gold", density=50)
}
}
