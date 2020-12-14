\name{inset}

\docType{methods}

\alias{inset,SpatVector-method}
\alias{inset,SpatRaster-method}

\title{Make an inset map}

\description{
Make an inset map
}

\usage{
\S4method{inset}{SpatVector}(x, e, loc="", scale=0.2, background="white", border="black", box=NULL, pbx, ...)

\S4method{inset}{SpatRaster}(x, e, loc="", scale=0.2, background="white", border="black", box=NULL, pbx, ...)
}

\arguments{
  \item{x}{SpatVector, SpatRaster}
  \item{e}{SpatExtent to set the size and location of the inset. Or missing}
  \item{loc}{character. One of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"}
  \item{scale}{numeric. The relative size of the inset, used when x is missing}
  \item{background}{color for the background of the inset. Use \code{NA} for no background color}
  \item{border}{color for the border around the inset. Use \code{NA} for no border}
  \item{box}{SpatExtent or missing, to draw a box on the inset}
  \item{pbx}{list with graphical arguments for the box}
  \item{...}{additional arguments passed to plot for the drawing of \code{x}}
}


\seealso{\code{\link{sbar}}, \code{\link{rescale}}, \code{\link{shift}}}

\value{scaled and shifted SpatVector or SpatRaster (invisibly)}

\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- v[v$NAME_2 == "Diekirch", ]

plot(x, density=10, col="blue")
inset(v)

e <- ext(c(6.2, 6.3, 49.9, 50))
b <- ext(x)+0.02
inset(v, e=e, col=cols, box=b)

cols <- rep("light grey", 12)
cols[2] <- "red"

inset(v, loc="bottomleft", col=cols, box=b, background="light yellow",
pbx=list(lwd=2, lty=2, col="blue"))

s <- inset(v, col=cols, loc="bottomright", scale=0.1)
s

}

\keyword{methods}
\keyword{spatial}

