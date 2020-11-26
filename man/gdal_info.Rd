\name{describe}

\alias{describe}
\alias{meta}
\alias{meta_sds}
\alias{meta,character-method}
\alias{describe_sds}
\alias{gdal_version}
\alias{gdal_warnings}

\title{describe}

\description{
Describe the properties of a file with spatial data (\code{describe}, \code{meta}) or the properties of the sub-datasets in file (\code{describe_sds} and \code{meta_sds}). 
}


\usage{
describe(x, options="", print=TRUE, open_opt="", ...)

describe_sds(x, ...)

\S4method{meta}{character}(x, ...)

meta_sds(x, parse=FALSE, ...)

gdal_version()

gdal_warnings(level)
}

\arguments{
  \item{x}{character (a filename)}
  \item{options}{character. A vector of valid options including "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct", "nofl", "checksum", "proj4", "listmdd", "mdd <value>" where <value> specifies a domain or 'all', "wkt_format <value>" where value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd <subdataset>" where <subdataset> is the name of a sub-dataset. See \url{https://gdal.org/programs/gdalinfo.html}}  
  \item{print}{logical. If \code{TRUE}, print the results}
  \item{parse}{logical. If \code{TRUE}, metadata is parsed into components}
  \item{open_opt}{character. Driver specific open options}   
  \item{level}{integer between 0 and 3 representing the type warnings from GDAL that is shown. 0 = warnings and errors; 1 = errors only; 2 = only irrecoverable errors; 3 = no warnings. The default setting is 2}
  \item{...}{additional arguments. None implemented}   

}


\value{
character (invisibly, if \code{print=FALSE})
}



\examples{
f <- system.file("ex/test.tif", package="terra")
describe(f)

#g <- describe(f, c("json", "nomd", "proj4"), print=FALSE)
#cat(g, "\n")
}

\keyword{spatial}

