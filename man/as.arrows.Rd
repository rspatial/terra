\name{as.arrows}

\alias{as.arrows}

\title{Add directional arrows from a SpatRaster to a plot}

\description{
Create directional arrows and add these to a plot, typically used to represent vector fields such as wind directions, drainage flow directions, or other angular data.
}

\usage{
as.arrows(x, unit=c("degrees", "radians", "flowdir"),
       arrow_length=NA,clockwise=FALSE, angle_from_h_axis=0,
       angle_from_x_axis=angle_from_h_axis, length=0,...)
}

\arguments{
  \item{x}{SpatRaster with values representing angles or directions}
  \item{unit}{Character. Unit of the angle values in \code{x}. One of \code{"degrees"}, \code{"radians"} or \code{"flowdir"}. See \code{\link{terrain}} for \code{flowdir}}
  \item{arrow_length}{Numeric. Length of the arrows to be plotted}
  \item{clockwise}{Logical. Default is \code{FALSE}. If \code{TRUE}, angles measured in radians or degrees are considered clockwise}
  \item{angle_from_x_axis}{Anti-clockwise angle from x (or horizontal) axis to the reference direction}
  \item{angle_from_h_axis}{Anti-clockwise angle from x (or horizontal) axis to the reference direction. In the case of northing wind direction it is often used \code{angle_from_x_axis=90, unit="deg", clockwise=TRUE}}
  \item{length,...}{Additional graphical parameters passed to \code{\link[graphics]{arrows}}}
}


\author{
Emanuele Cordano
}

\examples{

r <- rast(system.file("ex/elev_vinschgau.tif", package="terra"))
r <- aggregate(r, fact=2,fun=min)
d <- terrain(r, "flowdir")

plot(r,col=terrain.colors(10))
as.arrows(d, unit="flowdir", col="black")
}

\seealso{
\code{\link[graphics]{arrows}}, \code{\link[terra]{terrain}}
}
