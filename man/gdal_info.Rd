\name{describe}

\alias{describe}
\alias{gdal_version}

\title{describe}

\description{
Describe the properties of a file with raster data. (using the "GDALinfo" tool)
}


\usage{
describe(filename, options="", print=TRUE, open_opt="", ...)
gdal_version()
}

\arguments{
  \item{filename}{character}
  \item{options}{character. A vector of valid options including "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct", "nofl", "checksum", "proj4", "listmdd", "mdd <value>" where <value> specifies a domain or 'all', "wkt_format <value>" where value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd <subdataset>" where <subdataset> is the name of a sub-dataset. See \url{https://gdal.org/programs/gdalinfo.html}}  
  \item{print}{logical. If \code{TRUE}, print the results}
  \item{open_opt}{character. Driver specific open options}   
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

