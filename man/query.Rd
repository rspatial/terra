\name{query}

\docType{methods}

\alias{query}
\alias{query,SpatVectorProxy-method}


\title{Query a SpatVectorProxy object}

\description{
Query a SpatVectorProxy to extract a subset
} 

\usage{
\S4method{query}{SpatVectorProxy}(x, start=1, n=nrow(x), vars=NULL, where=NULL, 
    extent=NULL, filter=NULL, sql=NULL, what="")
}

\arguments{
\item{x}{SpatVectorProxy}
\item{start}{positive integer. The record to start reading at}
\item{n}{positive integer. The number of records requested}
\item{vars}{character. Variable names. Must be a subset of \code{names(x)}}
\item{where}{character. expression like "NAME_1='California' AND ID > 3" , to subset records. Note that start and n are applied after executing the where statement}
\item{extent}{Spat* object. The extent of the object is used as a spatial filter to select the geometries to read. Ignored if \code{filter} is not \code{NULL}}
\item{filter}{SpatVector. Used as a spatial filter to select geometries to read (the convex hull is used for lines or points)}
\item{sql}{character. Arbitrary SQL statement. If used, arguments "start", "n", "vars" and "where" are ignored}
\item{what}{character indicating what to read. Either \code{""} for geometries and attributes, or \code{"geoms"} to only read the geometries, \code{"attributes"} to only read the attributes (that are returned as a data.frame)}
}

\seealso{\code{\link{vect}}}


\value{SpatVector}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f, proxy=TRUE)
v

x <- query(v, vars=c("ID_2", "NAME_2"), start=5, n=2)
x

query(v, vars=c("ID_2", "NAME_1", "NAME_2"), where="NAME_1='Grevenmacher' AND ID_2 > 6")

## with an extent
e <- ext(5.9, 6.3, 49.9, 50) 
x <- query(v, extent=e)

## with polygons
p <- as.polygons(e)
x <- query(v, filter=p)
x
}

\keyword{methods}
\keyword{spatial}

