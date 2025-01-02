\name{extractRange}

\docType{methods}

\alias{extractRange}
\alias{extractRange,SpatRaster-method}
\alias{extractRange,SpatRaster,ANY-method}

\title{Extract values for a range of layers from a SpatRaster}

\description{
Extract values from a SpatRaster for a set of locations and a range of layers. To extract values for a single or all layers, use \code{\link{extract}}
} 

\usage{
\S4method{extractRange}{SpatRaster}(x, y, first, last, lyr_fun=NULL, 
		geom_fun=NULL, ID=FALSE, na.rm=TRUE, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{y}{SpatVector (points, lines, or polygons). Alternatively, for points, a 2-column matrix or data.frame (x, y) or (lon, lat). Or a vector with cell numbers}
 \item{first}{layer name of number, indicating the first layer in the range of layers to be considered}
\item{last}{layer name or number, indicating the last layer in the range to be considered}  
\item{lyr_fun}{function to summarize the extracted data across layers}
\item{geom_fun}{function to summarize the extracted data for each line or polygon geometry. Ignored if \code{y} has point geometry}
\item{ID}{logical. Should an ID column be added? If so, the first column returned has the IDs (record numbers) of \code{y}}
\item{na.rm}{logical. Should missing values be ignored?}
\item{...}{additional arguments passed to \code{extract}}
}

\value{numeric or data.frame}

\seealso{\code{\link{extract}}}

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
xy <- data.frame(c(50,80), c(30, 60))
extract(r, xy)
extract(r, xy, layer=c("red", "green"))
extractRange(r, xy, first=1:2, last=3:2, lyr_fun=sum)
}

\keyword{methods}
\keyword{spatial}

