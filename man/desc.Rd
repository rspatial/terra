\name{describe}

\alias{desc}
\alias{meta}
\alias{desc,character-method}
\alias{meta,character-method}

\title{describe}

\description{
Describe the properties of spatial data in a file, or get the metadata. The output from \code{desc} is generated with "GDALinfo".
}


\usage{
\S4method{desc}{character}(x, sds=FALSE, options="", print=FALSE, open_opt="", ...)

\S4method{meta}{character}(x, sds=FALSE, parse=FALSE, ...)
}

\arguments{
  \item{x}{character. The name of a file with spatial data. Or a fully specified subdataset within a file such as \code{"NETCDF:\"AVHRR.nc\":NDVI"}}
  \item{sds}{logical. If \code{TRUE} the description or metadata of the subdatasets is returned (if available)}
  \item{options}{character. A vector of valid options including "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct", "nofl", "checksum", "proj4", "listmdd", "mdd <value>" where <value> specifies a domain or 'all', "wkt_format <value>" where value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd <subdataset>" where <subdataset> is the name or identifier of a sub-dataset. See \url{https://gdal.org/programs/gdalinfo.html}. Ignored if \code{sds=TRUE}}  
  \item{print}{logical. If \code{TRUE}, print the results}
  \item{open_opt}{character. Driver specific open options}   
  \item{parse}{logical. If \code{TRUE}, metadata for subdatasets is parsed into components}
  \item{...}{additional arguments. None implemented}   
}


\value{
character (invisibly, if \code{print=FALSE})
}



\examples{
f <- system.file("ex/test.tif", package="terra")
desc(f)
meta(f)
#g <- desc(f, options=c("json", "nomd", "proj4"))
#cat(g, "\n")
}

\keyword{spatial}

