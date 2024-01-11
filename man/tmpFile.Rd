\name{tmpFiles}

\alias{tmpFiles}

\title{Temporary files}

\description{
List and optionally remove temporary files created by the terra package. These files are created when an output SpatRaster may be too large to store in memory (RAM). This can happen when no filename is provided to a function and when using functions where you cannot provide a filename.

Temporary files are automatically removed at the end of each R session that ends normally. You can use \code{tmpFiles} to see the files in the current sessions, including those that are orphaned (not connect to a SpatRaster object any more) and from other (perhaps old) sessions, and remove all the temporary files. 
}

\usage{
tmpFiles(current=TRUE, orphan=FALSE, old=FALSE, remove=FALSE)
}

\arguments{
  \item{current}{logical. If \code{TRUE}, temporary files from the current R session are included}
  \item{orphan}{logical. If \code{TRUE}, temporary files from the current R session that are no longer associated with a SpatRaster (if \code{current} is \code{TRUE} these are also included)}
  \item{old}{logical. If \code{TRUE}, temporary files from other "R" sessions. Unless you are running multiple instances of R at the same time, these are from old (possibly crashed) R  sessions and should be removed}
  \item{remove}{logical. If \code{TRUE}, temporary files are removed}
}


\value{
character
}


\seealso{ \code{\link{terraOptions}} }  

\examples{
tmpFiles()
}

\keyword{ spatial }


