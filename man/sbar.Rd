\name{sbar}

\alias{sbar}
\alias{arrow}

  
\title{scale bar and North arrow}

\description{
Add a scale bar or a North arrow to a map
}

\usage{
sbar(d, xy=NULL, type="line", divs=2, below="", 
    lonlat=NULL, label, adj=c(0.5, -1), lwd=2, xpd=TRUE, ...)

arrow(d, xy=NULL, type=1, label="N", angle=0, head=0.1, xpd=TRUE, ...)
}

\arguments{
\item{d}{numeric. Distance covered by the scale bar or arrow. For the scale bar, it should be in the units of the coordinates of the plot (map), and in km for angular (longitude/latitude) data; see argument \code{lonlat}. For the arrow it should be in plot coordinates}
\item{xy}{numeric. x and y coordinate to place the scale bar or arrow. It can also be one of following character values: "bottomleft", "bottom", "bottomright", topleft", "top", "topright", "left", "right", or NULL}
\item{type}{for \code{sbar} character: "line" or "bar". For \code{arrow}: integer between 1 and 12, or a character (unicode) representation of a right pointing arrow such as \code{"\u27A9"}}
\item{divs}{number of divisions for a bar: 2 or 4}
\item{below}{character. Text to go below the scale bar (e.g., "kilometers")}
\item{lonlat}{logical or NULL. If logical, \code{TRUE} indicates if the plot is using longitude/latitude coordinates. If \code{NULL} this is guessed from the plot's coordinates}
\item{adj}{adjustment for text placement}
\item{label}{vector of three numbers to label the scale bar (beginning, midpoint, end); or a character printed over the arrow}
\item{lwd}{line width for the "line" type of the scale bar}
\item{xpd}{logical. If \code{TRUE}, the scale bar or arrow can be outside the plot area}

\item{angle}{numeric. The angle of the arrow in degrees}
\item{head}{numeric. The size of the arrow "head", for \code{type=1}}
\item{...}{graphical arguments to be passed to other methods }
}

\value{
none
}

\seealso{ \code{\link[terra]{plot}}, \code{\link[terra]{inset}} }

\examples{
f <- system.file("ex/meuse.tif", package="terra")
r <- rast(f)
plot(r)
sbar()
arrow()
sbar(1000, xy=c(178500, 333500), type="bar", divs=4, cex=.8)
sbar(1000, xy="bottomright", divs=4, cex=.8)
arrow(250, c(178550, 332500))


f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
plot(r, type="interval")
sbar(20, c(6.2, 50.1), type="bar", cex=.8, divs=4)
sbar(15, c(6.3, 50), type="bar", below="km", label=c(0,7.5,15), cex=.8)
sbar(15, c(6.65, 49.8), cex=.8, label=c(0,"km",15))
arrow()
arrow(d=.05, xy=c(6.7, 49.97), angle=45, label="NE")

sbar(15, c(6.65, 49.7), cex=.8, label="15 kilometer", lwd=5)
sbar(15, c(6.65, 49.6), divs=4, cex=.8, below="km")

arrow(d=.05, xy=c(6.7, 49.5), angle=180, label="S", lwd=2, col="blue")

}


\keyword{spatial}
