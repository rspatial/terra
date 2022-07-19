\name{plet}

\docType{methods}

\alias{plet}
\alias{plet,SpatRaster-method}
\alias{plet,SpatVector-method}
\alias{lines,leaflet-method}
\alias{points,leaflet-method}

\title{Plot with leaflet}

\description{
Plot the values of a SpatRaster or SpatVector to make an interactive leaflet map that is displayed in a browser.

These methods require that \code{packageVersion("leaflet") > "2.1.1"} or the development version of leaflet that you can install with \code{remotes::install_github("rstudio/leaflet")}.
}

\usage{
\S4method{plet}{SpatRaster}(x, y=1, col, alpha=0.8, tiles=NULL, maxcell=500000, legend="bottomright", shared=FALSE, main=y, collapse=TRUE, map=NULL)

\S4method{plet}{SpatVector}(x, y="", col, alpha=1, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), legend="bottomright", main=names(x), popup=TRUE, split=FALSE, collapse=FALSE, cex=1, map=NULL)

\S4method{lines}{leaflet}(x, y, col, lwd=3, alpha=1)

\S4method{points}{leaflet}(x, y, col, cex=1, alpha=1, popup=FALSE)
}

\arguments{
  \item{x}{SpatRaster, SpatVector, or leaflet object}
  \item{y}{missing, or positive integer, or character (variable or layer name) indicating the layer(s) to be plotted. If \code{x} is a SpatRater, you can select multiple layers}
  \item{col}{character. Vector of colors or color generating function}
  \item{alpha}{Number between 0 and 1 to set the transparency for all colors (0 is transparent, 1 is opaque)}
  \item{tiles}{character or NULL. Names of background tile providers}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}  
  \item{legend}{character to indicate the legend position ("bottomleft", "bottomright", "topleft" or "topright") or NULL to suppress the legend}
  \item{main}{character. Title for the legend. The length should be 1 if \code{x} is a SpatVector and length nlyr(x) if \code{x} is a SpatVector}
  
  \item{shared}{logical. Should the legend be the same for all rasters (if multiple layers of SpatRaster \code{x} are mapped)}  
  \item{map}{leaflet object}
  \item{collapse}{logical. Should the layers "control" panel be collapsed?}
  \item{split}{logical. IF \code{TRUE} a check-box is created to toggle each value in code{y} (If \code{x} is a SpatVector)}
  \item{cex}{numeric. point size magnifier. See \code{\link{par}}}
  \item{lwd}{numeric, line-width. See \code{\link{par}}}
  \item{popup}{logical. Should pop-ups be created?}
}

\seealso{
 \code{\link{plot}}
}

\examples{ 
\dontrun{
if (require(leaflet)) {
if (packageVersion("leaflet") > "2.1.1") {

v <- vect(system.file("ex/lux.shp", package="terra"))
p <- spatSample(as.polygons(v, ext=T), 10)
values(p) = data.frame(id=11:20, name=letters[1:10])

m <- plet(v, "NAME_1", split=TRUE, alpha=.2)
m <- points(m, p, col="gray", cex=2, popup=T)
lines(m, v)

r <- rast(system.file("ex/elev.tif", package="terra"))
plet(r) |> lines(v, lwd=1)

plet(r, tiles="Streets") |> lines(v, lwd=2, col="blue")

x <- c(r, 50*classify(r, 5))
names(x) <- c("first", "second")

plet(x, 1:2, tiles="Streets", collapse=FALSE) |> lines(v, lwd=2, col="blue")
plet(x, 1:2, tiles="Streets", shared=TRUE, collapse=FALSE) |> lines(v, lwd=2, col="blue")

}}}}

\keyword{methods}
\keyword{spatial}

