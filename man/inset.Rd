\name{inset}

\docType{methods}

\alias{inset}
\alias{inset,SpatVector-method}
\alias{inset,SpatRaster-method}

\alias{inext}
\alias{inext,SpatVector-method}

\title{Make an inset map}

\description{
Make an inset map or scale the extent of a SpatVector
}

\usage{
\S4method{inset}{SpatVector}(x, e, loc="", scale=0.2, background="white", 
	perimeter=TRUE, box=NULL, pper, pbox, offset=0.1, add=TRUE, ...)

\S4method{inset}{SpatRaster}(x, e, loc="", scale=0.2, background="white", 
	perimeter=TRUE, box=NULL, pper, pbox, offset=0.1, add=TRUE, ...)

\S4method{inext}{SpatVector}(x, e, y=NULL, gap=0)
}

\arguments{
  \item{x}{SpatVector, SpatRaster}
  \item{e}{SpatExtent to set the size and location of the inset. Or missing}
  \item{loc}{character. One of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"}
  \item{scale}{numeric. The relative size of the inset, used when x is missing}
  \item{background}{color for the background of the inset. Use \code{NA} for no background color}
  \item{perimeter}{logical. If \code{TRUE} a perimeter (border) is drawn around the inset}
  \item{box}{SpatExtent or missing, to draw a box on the inset, e.g. to show where the map is located in a larger area}
  \item{pper}{list with graphical parameters (arguments) such as \code{col} and \code{lwd} for the perimeter line}
  \item{pbox}{list with graphical parameters (arguments) such as \code{col} and \code{lwd} for the box (line)}
  \item{offset}{numeric. Value between 0.1 and 1 to indicate the relative distance between what is mapped and the bounding box}
  \item{add}{logical. Add the inset to the map?}
  \item{...}{additional arguments passed to plot for the drawing of \code{x}}
  \item{y}{SpatVector. If not NULL, \code{y} is scaled based with the parameters for \code{x}. This is useful, for example, when \code{x} represent boundaries, and \code{y} points within these boundaries}
  \item{gap}{numeric to add space between the SpatVector and the SpatExtent}
}


\seealso{\code{\link{sbar}}, \code{\link{rescale}}, \code{\link{shift}}}


\value{scaled and shifted SpatVector or SpatRaster (returned invisibly)}

\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- v[v$NAME_2 == "Diekirch", ]

plot(x, density=10, col="blue")
inset(v)

# more elaborate
plot(x, density=10, col="blue")
inset(v, col = "brown", border="lightgrey", perimeter=TRUE, 
		pper=list(col="orange", lwd=3, lty=2), 
		box=ext(x), pbox=list(col="blue", lwd=2))

cols <- rep("light grey", 12)
cols[2] <- "red"
e <- ext(c(6.2, 6.3, 49.9, 50))
b <- ext(x)+0.02
inset(v, e=e, col=cols, box=b)

# with a SpatRaster
ff <- system.file("ex/elev.tif", package="terra")
r <- rast(ff)
r <- crop(r, ext(x) + .01)
plot(r, type="int", mar=c(2,2,2,2), plg=list(x="topright"))
lines(v, lwd=1.5)
lines(x, lwd=2.5)
inset(v, col=cols, loc="topleft", scale=0.15)

# a more complex one
plot(r, plg=list(title="meter\n", shrink=.2, cex=.8))
lines(v, lwd=4, col="white")
lines(v, lwd=1.5)
lines(x, lwd=2.5)
text(x, "NAME_2", cex=1.5, halo=TRUE)
sbar(6, c(6.04, 49.785), type="bar", below="km", label=c(0,3,6), cex=.8)
s <- inset(v, col=cols, box=b, scale=.2, loc="topright", background="light yellow",
	pbox=list(lwd=2, lty=5, col="blue"))

# note the returned inset SpatVector
s
lines(s, col="orange")
}

\keyword{methods}
\keyword{spatial}

