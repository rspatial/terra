\name{crosstab}

\docType{methods}

\alias{crosstab}
\alias{crosstab,SpatRaster,missing-method}

\title{Cross-tabulate}

\description{
Cross-tabulate the layers of a SpatRaster to create a contingency table.
}

\usage{
\S4method{crosstab}{SpatRaster,missing}(x, digits=0, long=FALSE, useNA=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{digits}{integer. The number of digits for rounding the values before cross-tabulation}
  \item{long}{logical. If \code{TRUE} the results are returned in 'long' format data.frame instead of a table}
  \item{useNA}{logical, indicting if the table should includes counts of \code{NA} values}
}


\value{
A table or data.frame
}

\seealso{ \code{\link{freq}}, \code{\link{zonal}} }



\examples{
r <- s <- rast(nc=5, nr=5)
set.seed(1)
values(r) <- runif(ncell(r)) * 2
values(s) <- runif(ncell(r)) * 3
x <- c(r, s)

crosstab(x)

rs <- r/s
r[1:5] <- NA
s[20:25] <- NA
x <- c(r, s, rs)
crosstab(x, useNA=TRUE, long=TRUE)
}

\keyword{methods}
\keyword{spatial}

