\name{click}

\alias{click}
\alias{click,SpatRaster-method}
\alias{click,missing-method}

\title{Query by clicking on a map}

\description{
Click on a map (plot) to get values of a SpatRaster at that location; and optionally the coordinates and cell number of the location. For SpatVectors you need to click twice (draw a box).
}


\usage{
\S4method{click}{SpatRaster}(x, n=Inf, id=FALSE, xy=FALSE, cell=FALSE, type="n", show=TRUE, ...)
}


\arguments{
  \item{x}{SpatRaster or SpatVector, or missing}
  \item{n}{number of clicks on the plot (map)}
  \item{id}{Logical. If \code{TRUE}, a numeric ID is shown on the map that corresponds to the row number of the output}
  \item{xy}{Logical. If \code{TRUE}, xy coordinates are included in the output}
  \item{cell}{Logical. If \code{TRUE}, cell numbers are included in the output}
  \item{type}{One of "n", "p", "l" or "o". If "p" or "o" the points are plotted; if "l" or "o" they are joined by lines. See ?locator}
  \item{show}{logical. Print the values after each click?}
  \item{...}{additional graphics parameters used if type != "n" for plotting the locations. See ?locator}   
 }
 
\value{
The value(s) of \code{x} at the point(s) clicked on (or touched by the box drawn).
} 

\note{
The plot only provides the coordinates for a spatial query, the values are read from the SpatRaster that is passed as an argument. Thus you can extract values from an object that has not been plotted, as long as it spatialy overlaps with with the extent of the plot.
  
Unless the process is terminated prematurely values at at most \code{n} positions are determined. The identification process can be terminated by clicking the second mouse button and selecting 'Stop' from the menu, or from the 'Stop' menu on the graphics window. 
}


\seealso{\link{draw}}


\examples{
\dontrun{
 r <- rast(system.file("external/test.grd", package="raster"))
 plot(r)
 click(r)
 # now click on the plot (map)
}}

\keyword{ spatial }
