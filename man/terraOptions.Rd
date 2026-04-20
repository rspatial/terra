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

\bold{datatype} - default data type. See \code{\link{writeRaster}}.

\bold{todisk} - logical. If \code{TRUE} write all raster data to disk (temp file if no file name is specified). For debugging.

\bold{progress} - non-negative integer. A progress bar is shown if the number of chunks in which the data is processed is larger than this number. No progress bar is shown if the value is zero.

\bold{verbose} - logical. If \code{TRUE} debugging info is printed for some functions.

\bold{tolerance} - numeric. Difference in raster extent (expressed as the fraction of the raster resolution) that can be ignored when comparing alignment of rasters.

\bold{parallel} - logical. The master switch for terra's internal parallelism. When \code{TRUE} (the default in version >= 1.9-21), the C++ kernels that have been parallelized with Intel TBB (parts of \code{arith}, distance calculations on \code{SpatVector}, the fast \code{focal} path, etc.) will use multiple threads. When \code{FALSE}, all such kernels run on a single thread. Use \code{terra:::.have_TBB()} to check whether TBB support was compiled in; without TBB, this option has no effect. The same option can be set per call as a write option, e.g. \code{focal(r, w=21, fun="mean", wopt=list(parallel=TRUE))}.

\bold{threads} - non-negative integer. Cap on the number of threads used by parallel kernels (when \code{parallel=TRUE}) and by GDAL warp (in \code{project}, \code{resample}, \code{warp}). \code{0} (the default) means "no cap": TBB picks a sensible default (usually all logical CPUs) and GDAL uses \code{NUM_THREADS=ALL_CPUS}. A positive value caps both at that count, which is useful to leave room for other processes or to keep forked workers from each spawning hundreds of children. The value is honoured by every TBB-parallel kernel in terra through a single \code{tbb::task_arena}, and is forwarded to GDAL warp's \code{NUM_THREADS} option.
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
terraOptions(parallel=TRUE, threads=4)   # use TBB, capped at 4 threads
terraOptions()
}

\keyword{classes}
\keyword{spatial}

