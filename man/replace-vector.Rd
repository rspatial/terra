\name{vector-attributes}

\alias{$}
\alias{$<-}

\alias{$,SpatVector-method}
\alias{$<-,SpatVector-method}

\title{Get or replace attribute values of a SpatVector}

\description{
Replace values of a SpatVector.  
}

\usage{
x$name
x$name <- value
}

\arguments{
  \item{x}{SpatVector}
  \item{name}{a literal character string (not a variable)}
  \item{value}{a vector}
}

\value{
vector 
}


\seealso{
\code{\link{values}}
}


\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
v$NAME_1
v$ID_1 <- LETTERS[1:12]
v$new <- sample(12)
values(v)
}

\keyword{ spatial }
