\name{subset-vector}

\alias{subset,SpatVector-method}

\alias{[,SpatVector,missing,missing-method}
\alias{[,SpatVector,missing,numeric-method}
\alias{[,SpatVector,missing,character-method}
\alias{[,SpatVector,missing,logical-method}
\alias{[,SpatVector,numeric,missing-method}
\alias{[,SpatVector,numeric,numeric-method}
\alias{[,SpatVector,numeric,character-method}
\alias{[,SpatVector,numeric,logical-method}
\alias{[,SpatVector,logical,missing-method}
\alias{[,SpatVector,logical,character-method}
\alias{[,SpatVector,logical,numeric-method}
\alias{[,SpatVector,logical,logical-method}
\alias{[,SpatVector,data.frame,missing-method}
\alias{[,SpatVector,matrix,missing-method}

\title{Subset of a SpatVector}

\description{
Select a subset of variables or records from a SpatVector.  
}

\usage{
\S4method{subset}{SpatVector}(x, subset, select, drop=FALSE)
}

\arguments{
  \item{x}{SpatVector}
  \item{subset}{logical expression indicating the rows to keep: missing values are taken as FALSE}
  \item{select}{expression, indicating columns to select} 
  \item{drop}{logical. If \code{TRUE}, the geometries will be dropped, and a data.frame is returned}
}

\value{
SpatVector or, if \code{drop=TRUE}, a data.frame.
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v[2:3,]
v[,2:3]
subset(v, v$NAME_1 == "Diekirch") 
}

\keyword{ spatial }
