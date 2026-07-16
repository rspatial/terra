\name{arrows_on_rast}
\alias{arrows_on_rast}
\title{Visualize Directional Arrows on a SpatRaster}
\description{
This function overlays directional arrows on a plotted \code{\link[terra]{SpatRaster}} object, 
typically used to represent vector fields such as wind directions, drainage flow directions, or other angular data.
}
\usage{
arrows_on_rast(x="flowdir",x_base=NULL,unit=c("deg","rad","flowdir"),
                arrow_length=NA,clockwise=FALSE,
                angle_from_h_axis=0,
                angle_from_x_axis=angle_from_h_axis,
                length=0,...)
}
\arguments{
  \item{x}{A \code{\link[terra]{SpatRaster}} object representing angles or directions (in degrees or radians or flowdirs).}
  \item{x_base}{Optional. A \code{\link[terra]{SpatRaster}} object to use as a basemap if no plot is already present.}
  \item{unit}{Character. Unit of the angle values in \code{x}. Either \code{"degrees"} or \code{"radians"} or \code{"flowdir"}. See \code{\link[terra]{terrain}} for \code{flowdir} option.}
  \item{arrow_length}{Numeric. Length of the arrows to be plotted.}
  \item{clockwise}{Logical. Default is \code{FALSE}. If \code{TRUE}, angles measured in radians or degrees are considered clockwise.}
  \item{angle_from_x_axis}{Anti-clockwise angle from x (or horizontal) axis to the reference direction. Default is 0.}
  \item{angle_from_h_axis}{Anti-clockwise angle from x (or horizontal) axis to the reference direction. Default is 0, e.g., in the case of northing wind direction it is often used \code{angle_x_axis=90, unit="deg", clockwise=TRUE}.}
  \item{length,...}{Additional graphical parameters passed to \code{\link[graphics]{arrows}}.}
}
\value{
Invisibly returns the coordinates and vector components used for plotting.
}
\author{
Emanuele Cordano
}
\examples{

f <- system.file("ex/elev_vinschgau.tif", package="terra")
r <- rast(f)  |> aggregate(fact=2,fun=min)
d <- terrain(r, "flowdir")


plot(r,col=terrain.colors(10))
arrows_on_rast(d, unit="flowdir",col="black")



}
\seealso{
\code{\link[graphics]{arrows}}, \code{\link[terra]{SpatRaster}}, \code{\link[terra]{terrain}}
}
