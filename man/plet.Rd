\name{plet}

\docType{methods}

\alias{plet}
\alias{plet,SpatRaster-method}
\alias{plet,SpatVector-method}
\alias{lines,leaflet-method}
\alias{points,leaflet-method}

\title{Make a leaflet map}

\description{
Plot the values of a SpatRaster or SpatVector to make an interactive leaflet map that is displayed in a browser

These methods require that \code{packageVersion("leaflet") > "2.1.1"} or the development version of leaflet that you can install with \code{remotes::install_github("rstudio/leaflet")}
}

\usage{
\S4method{plet}{SpatRaster}(x, y=1, col=rainbow, split=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), alpha=1, legend="bottomright", collapse=FALSE, cex=1, map=NULL)

\S4method{plet}{SpatVector}(x, y="", col, alpha=0.2, tiles="", maxcell=500000, legend="bottomright", map=NULL)

\S4method{lines}{leaflet}(x, y, col, lwd=3, alpha=1)

\S4method{points}{leaflet}(x, y, col, cex=1, alpha=1, popup=FALSE)
}

\arguments{
  \item{x}{SpatRaster, SpatVector, or leaflet object}
  \item{y}{missing or positive integer or name indicating the layer(s) to be plotted}
  \item{col}{character. Colors. The default is \code{rev(grDevices::terrain.colors(50))}}
  \item{legend}{character to indicate the legend position ("bottomleft", "bottomright", "topleft" or "topright") or NULL to suppress the legend}
  \item{tiles}{character. Names of background tile providers}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plotlet}
  \item{alpha}{Number between 0 and 1 to set the transparency for all colors (0 is transparent, 1 is opaque)}
  \item{map}{leaflet object}
}

\seealso{
 \code{\link{plot}}
}

\examples{ 
#v <- vect(system.file("ex/lux.shp", package="terra"))
#p <- spatSample(as.polygons(v, ext=T), 10)
#values(p) = data.frame(id=11:20, name=letters[1:10])

#m <- plotlet(v, "NAME_1", split=TRUE, fillOpacity=.2)
#m <- points(m, p, col="gray", cex=2, popup=T)
lines(m, v)

#r <- rast(system.file("ex/elev.tif", package="terra"))
#plotlet(r) |> lines(v, lwd=1)

#plotlet(r, tiles="OpenWeatherMap.Wind") |> lines(v, lwd=2, col="blue")
}


\keyword{methods}
\keyword{spatial}

