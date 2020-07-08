\name{options}


\alias{terraOptions}


\title{Options}

\description{
Class and methods for showing and setting general options for terra.
}

\usage{
terraOptions(...)
}

\arguments{
  \item{...}{additional arguments. Can be use to set options (see examples). When empty, the current options are shown}  
} 

\details{
The following options are available. 

memfrac - value between 0.1 and 0.8. The fraction of RAM that may be used by the program.

tempdir - directory where temporary files are written. The default what is returned by tempdir().

datatype - default data type. See \code{\link{writeRaster}}
  
filetype - default file type. See \code{\link{writeRaster}}

progress - non-negative integer. A progress bar is shown if the number of chunks in which the data is processed is larger than this number. No progress bar is shown if the value is zero

verbose - logical. If \code{TRUE} debugging info is printed for some functions}
}

\examples{
terraOptions()
terraOptions(memfrac=0.5, tempdir = "c:/temp")
terraOptions(progress=10)
terraOptions()
}


\keyword{classes}
\keyword{spatial}
