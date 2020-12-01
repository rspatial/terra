\name{NAflag}

\alias{NAflag<-}
\alias{NAflag}
\alias{NAflag,SpatRaster-method}
\alias{NAflag<-,SpatRaster-method}


\title{Set a value to NA}

\description{
Set a particular value of a SpatRaster to \code{NA}. If the values are in memory the change is made immediately. If the values are on disk, the changes is only made for values read from the file ("lazy evaluation"). 
}

\usage{
\S4method{NAflag}{SpatRaster}(x, ...)

\S4method{NAflag}{SpatRaster}(x)<-value
}

\arguments{
 \item{x}{SpatRaster}
 \item{value}{numeric. The value to be interpreted as NA; set this before reading the values from the file. This can be a single value, or multiple values, one for each data source (file / subdataset)}
 \item{...}{additional arguments. None implemented}
}


\value{
none or numeric
}


\examples{ 
s <- rast(system.file("ex/logo.tif", package="terra"))[[1]]   
NAflag(s) <- 255
plot(s)
NAflag(s)
}

\keyword{ spatial }

