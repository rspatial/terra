\name{colors}

\docType{methods}

\alias{has.colors}
\alias{has.colors,SpatRaster-method}
\alias{coltab}
\alias{coltab,SpatRaster-method}
\alias{coltab<-}
\alias{coltab<-,SpatRaster-method}
  
\title{Color table}

\description{
Get or set color table(s) associated with a SpatRaster. Color tables are used for associating colors with values, for use in mapping (plot).
}

\usage{
\S4method{coltab}{SpatRaster}(x)

\S4method{coltab}{SpatRaster}(x, ..., layer=1)<-value

\S4method{has.colors}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{a two-column data.frame (first column the cell value, the second column the color); 
  a vector of colors (the first one is the color for value 0 and so on); or a four (value,red,green,blue) or five (including alpha) column data.frame also from 0 to n; or NULL to remove the color table. You can also supply a list of such data.frames to set a color table to all layers}
  \item{...}{additional arguments (none implemented)}
}  


\value{
data.frame
}

\examples{
r <- rast(ncols=3, nrows=2, vals=1:6)
coltb <- data.frame(value=1:6, col=rainbow(6, end=.9))
coltb

plot(r)

has.colors(r)
coltab(r) <- coltb
plot(r)
has.colors(r)

tb <- coltab(r)
class(tb)
dim(tb[[1]])
}

\keyword{methods}
\keyword{spatial}

