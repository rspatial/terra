\name{svc}

\docType{methods}

\alias{svc}
\alias{svc,missing-method}
\alias{svc,character-method}
\alias{svc,list-method}
\alias{svc,sf-method}
\alias{svc,SpatVector-method}
\alias{[<-,SpatVectorCollection,numeric,missing-method}
\alias{[,SpatVectorCollection,numeric,missing-method}
\alias{[[,SpatVectorCollection,numeric,missing-method}



\title{Create a SpatVectorCollection}

\description{
Methods to create a SpatVectorCollection. This is an object to hold "sub-datasets", each a SpatVector, perhaps of different geometry type. 
}

\usage{
\S4method{svc}{missing}(x) 

\S4method{svc}{SpatVector}(x, ...) 

\S4method{svc}{list}(x) 

\S4method{svc}{character}(x, layer="", query="", extent=NULL, filter=NULL)
}

\arguments{
\item{x}{SpatVector, character (filename), list with SpatVectors, or missing}
\item{...}{Additional \code{SpatVector}s} 
\item{layer}{character. layer name to select a layer from a file (database) with multiple layers}
\item{query}{character. An query to subset the dataset in the \href{https://gdal.org/user/ogr_sql_dialect.html}{OGR-SQL dialect}}
\item{extent}{Spat* object. The extent of the object is used as a spatial filter to select the geometries to read. Ignored if \code{filter} is not \code{NULL}}
\item{filter}{SpatVector. Used as a spatial filter to select geometries to read (the convex hull is used for lines or points). It is guaranteed that all features that overlap with the extent of filter will be returned. It can happen that additional geometries are returned}
}


\value{
SpatVectorCollection
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- svc()
x <- svc(v, v[1:3,], as.lines(v[3:5,]), as.points(v))
length(x)
x

# extract
x[3]

# replace
x[2] <- as.lines(v[1,])


}


\keyword{methods}
\keyword{spatial}

