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
  \item{...}{option names and values (see Details). Or missing, to show the current options}  
} 

\details{
The following options are available. 

\bold{memfrac} - value between 0 and 0.9 (larger values give a warning). The fraction of RAM that may be used by the program.

\bold{memmax} - the maximum amount of RAM (in GB) that terra can use when processing a raster dataset. Should be less than what is detected (see \code{\link{mem_info}}, and higher values are ignored. Set it to a negative number or NA to not set this option). \code{terraOptions} only shows the value of \code{memmax} it it set.

\bold{tempdir} - directory where temporary files are written. The default what is returned by \code{tempdir()}.

\bold{datatype} - default data type. See \code{\link{writeRaster}}

\bold{todisk} - logical. If \code{TRUE} write all raster data to disk (temp file if no file name is specified). For debugging.

\bold{progress} - non-negative integer. A progress bar is shown if the number of chunks in which the data is processed is larger than this number. No progress bar is shown if the value is zero

\bold{verbose} - logical. If \code{TRUE} debugging info is printed for some functions
}

\examples{
terraOptions()
terraOptions(memfrac=0.5, tempdir = "c:/temp")
terraOptions(progress=10)
terraOptions()
}

\keyword{classes}
\keyword{spatial}
