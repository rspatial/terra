\name{click}

\alias{click}
\alias{click,SpatRaster-method}
\alias{click,SpatVector-method}
\alias{click,missing-method}

\title{Query by clicking on a map}

\description{
Click on a map (plot) to get the coordinates or the values of a SpatRaster or SpatVector at that location. For a SpatRaster you can also get the coordinates and cell number of the location. 

This does to work well on the default RStudio plotting device. To work around that, you can first run \code{dev.new(noRStudioGD = TRUE)} which will create a separate window for plotting, then use \code{plot()} followed by \code{click()} and click on the map .
}


\usage{
\S4method{click}{SpatRaster}(x, n=10, id=FALSE, xy=FALSE, cell=FALSE, type="p", show=TRUE, ...)

\S4method{click}{SpatVector}(x,  n=10, id=FALSE, xy=FALSE, type="p", show=TRUE, ...)

\S4method{click}{missing}(x, n=10, id=FALSE, type="p", show=TRUE, ...)
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
A \code{data.frame} with the value(s) of all layers of SpatRaster \code{x} for the cell(s) clicked on; or with the attributes of the geometries of SpatVector \code{x} that intersect with the box drawn).
} 

\note{
The plot only provides the coordinates for a spatial query, the values are read from the SpatRaster or SpatVector that is passed as an argument. Thus you can extract values from an object that has not been plotted, as long as it spatially overlaps with with the extent of the plot.
  
Unless the process is terminated prematurely values at at most \code{n} positions are determined. The identification process can be terminated, depending on how you interact with R, by hitting Esc, or by clicking the right mouse button and selecting "Stop" from the menu, or from the "Stop" menu on the graphics window. 
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
