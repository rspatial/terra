\name{map.pal}

\alias{map.pal}
  
\title{color palettes for mapping}

\description{
Get a color palette for mapping. These palettes were copied from GRASS.
}

\usage{
map.pal(name, n=50, ...)

}

\arguments{
\item{name}{character (name of a palette, see Details), or missing (to get the available names)}
\item{n}{numeric. The number of colors}
\item{...}{additional arguments that are passed to \code{\link{colorRamp}}}
}

\value{
none
}

\details{

\tabular{ll}{
Name \tab Description \cr
aspect \tab aspect oriented grey colors\cr
bcyr \tab blue through cyan through yellow to red\cr
bgyr \tab blue through green through yellow to red\cr
blues \tab white to blue\cr
byg \tab blue through yellow to green\cr
byr \tab blue through yellow to red\cr
curvature \tab for terrain curvatures\cr
differences \tab differences oriented colors\cr
elevation \tab maps relative ranges of raster values to elevation color ramp\cr
grass \tab GRASS GIS green (perceptually uniform)\cr
greens \tab white to green\cr
grey \tab grey scale\cr
gyr \tab green through yellow to red\cr
haxby \tab relative colors for bathymetry or topography\cr
inferno \tab perceptually uniform sequential colors inferno\cr
magma \tab perceptually uniform sequential colors\cr
oranges \tab white to orange\cr
plasma \tab perceptually uniform sequential colors\cr
rainbow \tab rainbow colors\cr
ramp \tab color ramp\cr
random \tab random colors\cr
reds \tab white to red\cr
roygbiv \tab \cr
rstcurv \tab terrain curvature\cr
ryb \tab red through yellow to blue\cr
ryg \tab red through yellow to green\cr
sepia \tab yellowish-brown through to white\cr
viridis \tab perceptually uniform sequential colors\cr
water \tab water depth\cr
wave \tab color wave\cr
}
}

\seealso{ \code{\link{terrain.colors}} }

\examples{
map.pal("elevation", 10)

r <- rast(system.file("ex/elev.tif", package="terra"))
plot(r, col=map.pal("elevation"))

map.pal()
}


\keyword{spatial}
