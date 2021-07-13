\name{vector-attributes}

\alias{$}
\alias{$<-}

\alias{$,SpatVector-method}
\alias{$<-,SpatVector-method}

\alias{[[,SpatVector,numeric,missing-method}
\alias{[[,SpatVector,logical,missing-method}
\alias{[[,SpatVector,character,missing-method}
\alias{[[<-,SpatVector,character,missing-method}
\alias{[[<-,SpatVector,numeric,missing-method}

\alias{[<-,SpatVector,ANY,ANY-method}
\alias{[<-,SpatVector,ANY,missing-method}
\alias{[<-,SpatVector,missing,ANY-method}

\title{Get or replace attribute values of a SpatVector}

\description{
Replace values of a SpatVector.  
}

\usage{
\S4method{$}{SpatVector}(x, name)

\S4method{$}{SpatVector}(x, name)<-value
}

\arguments{
  \item{x}{SpatVector}
  \item{name}{character (field name) or numeric (column number}
  \item{value}{vector of new values}
}

\value{
vector 
}


\seealso{
\code{\link{values}}
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$NAME_1
v$NAME_1[3] <- "my name"
v$ID_1 <- LETTERS[1:12]
v$new <- sample(12)
values(v)

v[2,2] <- "hello"
v[1,] <- v[10,]
v[,3] <- v[,1]
v[2, "NAME_2"] <- "terra"
head(v, 3)
}

\keyword{ spatial }
