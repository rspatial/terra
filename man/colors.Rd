\name{colors}

\docType{methods}

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

\S4method{coltab}{SpatRaster}(x, layer=1)<-value
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{a two-column data.frame (first column the cell value, second the color); 
  a vector of colors (the first one is the color for value 0 and so on); or a three (red,green,blue) or four (alpha) column data.frame also from 0 to n; or NULL to remove the color table}
}  


\value{
data.frame
}

\examples{
r <- rast(ncols=3, nrows=2, vals=0:5)
coltb <- data.frame(t(col2rgb(rainbow(6, end=.9), alpha=TRUE)))
coltb

plot(r)
coltab(r) <- coltb
plot(r)

tb <- coltab(r)
class(tb)
dim(tb[[1]])
}

\keyword{methods}
\keyword{spatial}

