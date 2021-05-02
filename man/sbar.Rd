\name{sbar}

\alias{sbar}
  
\title{scalebar}

\description{
Add a scalebar to a plot
}

\usage{
sbar(d, xy=NULL, type="line", divs=2, below="", 
    lonlat=NULL, label, adj=c(0.5, -1), lwd=2, xpd=TRUE, ...)
}

\arguments{
\item{d}{numeric. Distance covered by the scalebar. In the units of the coordinates of the plot, and in km for angular (longitude/latitude) data; see \code{lonlat}}
\item{xy}{x and y coordinate to place the plot. Can be NULL. Use \code{xy=click()} to make this interactive }
\item{type}{"line" or "bar"}
\item{divs}{number of divisions for a bar: 2 or 4}
\item{below}{character. Text to go below scalebar (e.g., "kilometers")}
\item{lonlat}{logical or NULL. If logical, \code{TRUE} indicates if the plot is using longitude/latitude coordinates. If \code{NULL} this is guessed from the plot's coordinates}
\item{adj}{adjustment for text placement}
\item{label}{vector of three numbers to label the scale bar (beginning, midpoint, end)}
\item{lwd}{line width for the "line" type scalebar}
\item{xpd}{logical. If \code{TRUE}, the scalebar can be (partly) outside the plot area}
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
sbar(1000)
sbar(1000, xy=c(178500, 333500), type="bar", divs=4, cex=.8)

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
plot(r, type="interval")
sbar(20, c(6.2, 50.1), type="bar", cex=.8, divs=4)
sbar(15, c(6.3, 50), type="bar", below="km", label=c(0,7.5,15), cex=.8)
sbar(15, c(6.65, 49.9), cex=.8, label=c(0,"km",15))

sbar(15, c(6.65, 49.8), cex=.8, label="15 kilometer", lwd=5)

sbar(15, c(6.65, 49.7), divs=4, cex=.8, below="km")
}

\keyword{spatial}
