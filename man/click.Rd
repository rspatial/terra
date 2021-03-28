\name{click}

\alias{click}
\alias{click,SpatRaster-method}
\alias{click,SpatVector-method}
\alias{click,missing-method}

\title{Query by clicking on a map}

\description{
Click on a map (plot) to get the coordinates or the values of a SpatRaster or SpatVector at that location. For a SpatRaster you can also get the coordinates and cell number of the location. 
}


\usage{
\S4method{click}{SpatRaster}(x, n=1000, id=FALSE, xy=FALSE, cell=FALSE, type="p", show=TRUE, ...)

\S4method{click}{SpatVector}(x,  n=1, id=FALSE, type="p", ...)

\S4method{click}{missing}(x, n=1, id=FALSE, type="p", ...)
}


\arguments{
  \item{x}{SpatRaster or SpatVector, or missing}
  \item{n}{number of clicks on the plot (map)}
  \item{id}{logical. If \code{TRUE}, a numeric ID is shown on the map that corresponds to the row number of the output}
  \item{xy}{logical. If \code{TRUE}, xy coordinates are included in the output}
  \item{cell}{logical. If \code{TRUE}, cell numbers are included in the output}
  \item{type}{one of "n", "p", "l" or "o". If "p" or "o" the points are plotted; if "l" or "o" they are joined by lines. See ?locator}
  \item{show}{logical. Print the values after each click?}
  \item{...}{additional graphics parameters used if type != "n" for plotting the locations. See ?locator}   
 }
 
\value{
The value(s) of \code{x} at the point(s) clicked on (or touched by the box drawn).
} 

\note{
The plot only provides the coordinates for a spatial query, the values are read from the SpatRaster that is passed as an argument. Thus you can extract values from an object that has not been plotted, as long as it spatially overlaps with with the extent of the plot.
  
Unless the process is terminated prematurely values at at most \code{n} positions are determined. The identification process can be terminated by hitting Esc, or by clicking the right mouse button and selecting "Stop" from the menu, or from the "Stop" menu on the graphics window. 
}


\seealso{\link{draw}}


\examples{
\dontrun{
r <-rast(system.file("ex/elev.tif", package="terra"))
plot(r)
click(r, n=1)
## now click on the plot (map)
}}

\keyword{ spatial }
