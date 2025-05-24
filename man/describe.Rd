\name{describe}

\alias{describe}
\alias{describe,character-method}
\alias{describe,SpatRaster-method}

\title{describe}

\description{
Describe the properties of spatial data in a file as generated with the "GDALinfo" tool.
}

\usage{
\S4method{describe}{character}(x, sds=FALSE, meta=FALSE, parse=FALSE, options="", print=FALSE, open_opt="", mdim=FALSE)

\S4method{describe}{SpatRaster}(x, source, ...)
}

\arguments{
  \item{x}{character. The name of a file with spatial data. Or a fully specified subdataset within a file such as \code{"NETCDF:\"AVHRR.nc\":NDVI"}}
  \item{sds}{logical. If \code{TRUE} the description or metadata of the subdatasets is returned (if available)}
  \item{meta}{logical. Get the file level metadata instead}
  \item{parse}{logical. If \code{TRUE}, metadata for subdatasets is parsed into components (if \code{meta=TRUE})}
  \item{options}{character. A vector of valid options (if \code{meta=FALSE}) including "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct", "nofl", "checksum", "proj4", "listmdd", "mdd <value>" where <value> specifies a domain or 'all', "wkt_format <value>" where value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd <subdataset>" where <subdataset> is the name or identifier of a sub-dataset. See \url{https://gdal.org/en/latest/programs/gdalinfo.html}. Ignored if \code{sds=TRUE}}  
  \item{print}{logical. If \code{TRUE}, print the results}
  \item{open_opt}{character. Driver specific open options}   
  \item{mdim}{logical. If \code{TRUE} return the GDAL multi-dimensional API description}}   
  \item{source}{positive integer between 1 and \code{nsrc(x)}}   
  \item{...}{additional arguments passed to the \code{describe<character>} method}
}

\value{
character (invisibly, if \code{print=FALSE})
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
describe(f)
describe(f, meta=TRUE)
#g <- describe(f, options=c("json", "nomd", "proj4"))
#head(g)
}

\keyword{spatial}
