\name{subset_dollar}

\alias{$}

\alias{$,SpatRaster-method}
\alias{$,SpatRasterDataset-method}
\alias{$,SpatVector-method}
\alias{$,SpatExtent-method}


\title{Subset a SpatRaster or a SpatVector}

\description{
Select a subset of layers from a SpatRaster or select a subset of records (row) and/or variables (columns) from a SpatVector.  
}

\usage{
\S4method{$}{SpatExtent}(x, name)
}

\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{name}{character. If \code{x} is a SpatRaster: layer name. If \code{x} is a SpatVector: variable name. If \code{x} is a SpatExtent: xmin, xmax, ymin or ymax}
}

\value{
if \code{x} is a \code{SpatRaster}: SpatRaster

if \code{x} is a \code{SpatVector}: SpatVector or, if \code{drop=TRUE}, a \code{data.frame}.
}


\seealso{\code{\link{subset}}, \code{\link{[}}, \code{\link{[[}}, \code{\link{extract}}}

\examples{
### SpatRaster
s <- rast(system.file("ex/logo.tif", package="terra"))   
subset(s, 2:3)
subset(s, c(3,2,3,1))
#equivalent to 
s[[ c(3,2,3,1) ]]

s[[c("red", "green")]]
s$red

# expression based (partial) matching of names with single brackets
s["re"]
s["^re"]

# not with double brackets
# s[["re"]]


### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[2:3,]
v[1:2, 2:3]

subset(v, v$NAME_1 == "Diekirch", c("NAME_1", "NAME_2"))

subset(v, NAME_1 == "Diekirch", c(NAME_1, NAME_2), NSE=TRUE)
}

\keyword{ spatial }
