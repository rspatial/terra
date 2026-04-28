\name{NAflag}

\alias{NAflag<-}
\alias{NAflag}
\alias{NAflag,SpatRaster-method}
\alias{NAflag<-,SpatRaster-method}


\title{Set the NA flag}

\description{
The main purpose of this method is to allow correct reading of a SpatRaster that is based on a file that has an incorrect or missing NA flag. The file is not changed, but the flagged value is set to NA when values are read from the file. In contrast, if the cell values are in memory the change is made immediately. 

To change cell values, you can use \code{\link{classify}}
}

\usage{
\S4method{NAflag}{SpatRaster}(x)

\S4method{NAflag}{SpatRaster}(x)<-value
}

\arguments{
 \item{x}{SpatRaster}
 \item{value}{numeric. The value to be interpreted as NA; set this before reading the values from the file. This can be a single value, or multiple values, one for each data source (file / subdataset)}
}


\value{
none or numeric
}

\seealso{
 \code{\link{classify}}
}

\examples{ 
s <- rast(system.file("ex/logo.tif", package="terra"))[[1]]   
NAflag(s) <- 255
plot(s)
NAflag(s)
}

\keyword{ spatial }

