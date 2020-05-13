\name{names}

\alias{names}
\alias{names<-}
\alias{names,SpatRaster-method}
\alias{names<-,SpatRaster-method}
\alias{names,SpatVector-method}
\alias{names<-,SpatVector-method}
  
\title{Names of Spat objects}

\description{
Get or set the names of the layers of a SpatRaster or the attributes of a SpatVector.
}

\usage{
\S4method{names}{SpatRaster}(x)

\S4method{names}{SpatRaster}(x)<-value

\S4method{names}{SpatVector}(x)

\S4method{names}{SpatVector}(x)<-value
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{value}{character (vector)}
}

\value{
character 
}


\examples{
s <- rast(ncols=5, nrows=5, nlyr=3)
nlyr(s)
names(s)
names(s) <- c("a", "b", "c")
names(s)

# space is not valid
names(s)[2] <- "hello world"
names(s)

# two invalid names
names(s) <- c("a", "  a  ", "3")
names(s)

# SpatVector names
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
names(v)
names(v) <- paste0(substr(names(v), 1, 2), "_", 1:ncol(v))
names(v)
}

\keyword{spatial}

