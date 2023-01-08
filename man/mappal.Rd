\name{map.pal}

\alias{map.pal}
  
\title{map.pal arrow}

\description{
Get a color palette for mapping. These palettes were copied from GRASS
}

\usage{
map.pal(name, n=50, ...)

}

\arguments{
\item{name}{character (name of a palette), or missing (to get the available names)}
\item{n}{numeric. The number of colors}
\item{...}{additional arguments that are passed to \code{\link{colorRamp}}}
}

\value{
none
}

\seealso{ \code{\link{terrain.colors}} }

\examples{
map.pal("elevation", 10)

r <- rast(system.file("ex/elev.tif", package="terra"))
plot(r, col=map.pal("elevation"))

map.pal()
}


\keyword{spatial}
