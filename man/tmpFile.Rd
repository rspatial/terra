\name{tmpFiles}

\alias{tmpFiles}
\alias{orphanTmpFiles}
\alias{staleTmpFiles}

\title{Temporary files}

\description{
List and optionally remove temporary files created by the terra package. These files are created when an output SpatRaster may be too large to store in memory (RAM). This can happen when no filename is provided to a function and in functions where you cannot provide a filename.

Temporary files are automatically removed at the end of each session that ends normally. You can use \code{tmpFiles} to see the files in the current sessions, from other (perhaps old) sessions, and remove all the temporary files. 
}

\usage{
tmpFiles(old=FALSE, remove=FALSE)
orphanTmpFiles()
staleTmpFiles()
}

\arguments{
  \item{old}{logical. If \code{TRUE}, temporary files from other (possibly aborted) sessions are included}
  \item{remove}{logical. If \code{TRUE}, temporary files are removed}
}


\value{
character
}


\seealso{ \code{\link{terraOptions}}, \code{\link[base]{tempfile}} }  

\examples{
tmpFiles()
}

\keyword{ spatial }


