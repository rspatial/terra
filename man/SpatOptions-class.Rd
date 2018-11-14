\name{SpatOptions-class}

\docType{class}

\alias{SpatOptions}
\alias{SpatOptions-class}
\alias{Rcpp_SpatOptions-class}
\alias{terraOptions}


\title{options}

\description{
Class and methods for showing and setting options.
}

\usage{
terraOptions(...)
}

\arguments{
  \item{...}{Additonal arguments. Use to set options. When empty, the current options are shown}  
}


\examples{
terraOptions()
terraOptions(memfrac=0.5, tempdir = "c:/temp")
terraOptions()
}


\keyword{classes}
\keyword{spatial}
