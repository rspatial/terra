\name{options}

\alias{terraOptions}

\title{Options}

\description{
Get or set general options.
}

\usage{
terraOptions(..., print=TRUE)
}

\arguments{
  \item{...}{option names and values (see Details). Or missing, to get or show the current options}  
  \item{print}{logical. If \code{TRUE} the option names and values are printed }  
} 

\details{
The following options are available. 

\bold{memfrac} - value between 0 and 0.9 (larger values give a warning). The fraction of RAM that may be used by the program.

\bold{memmin} - if memory required is below this threshold (in GB), the memory is assumed to be available. Otherwise, terra checks if it is available.

\bold{memmax} - the maximum amount of RAM (in GB) that terra is allowed to use when processing a raster dataset. Should be less than what is detected (see \code{\link{mem_info}}), and higher values are ignored. Set it to a negative number or NA to not set this option. \code{terraOptions} only shows the value of \code{memmax} if it is set.

\bold{tempdir} - directory where temporary files are written. The default what is returned by \code{tempdir()}.

\bold{datatype} - default data type. See \code{\link{writeRaster}}

\bold{todisk} - logical. If \code{TRUE} write all raster data to disk (temp file if no file name is specified). For debugging.

\bold{progress} - non-negative integer. A progress bar is shown if the number of chunks in which the data is processed is larger than this number. No progress bar is shown if the value is zero

\bold{verbose} - logical. If \code{TRUE} debugging info is printed for some functions
}

\note{
It is possible to set your own default options in "etc/.Rprofile.site" of your R installation like this

\code{options(terra_default=list(tempdir="d:/temp", memfrac=.4))}

But that may not be a good practice. It is clearer to set your favorite options at the beginning of each script.
}

\value{
list. Invisibly if \code{print=TRUE}
}

\examples{
terraOptions()
terraOptions(memfrac=0.5, tempdir = "c:/temp")
terraOptions(progress=10)
terraOptions()
}

\keyword{classes}
\keyword{spatial}

