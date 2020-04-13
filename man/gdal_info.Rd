\name{gdal_info}

\alias{gdal_info}

\title{GDALinfo}

\description{
This function provides a wrapper to the GDALinfo tool that provides the properties of a file with raster data. 
}


\usage{
gdal_info(filename, options="", print=TRUE, open_opt="", ...)
}

\arguments{
  \item{filename}{character}
  \item{options}{character. A vector of valid options including "json", "mm", "stats", "hist", "nogcp", "nomd", "norat", "noct", "nofl", "checksum", "proj4", "listmdd", "mdd <value>" where <value> specifies a domain or 'all', "wkt_format <value>" where value is one of 'WKT1', 'WKT2', 'WKT2_2015', or 'WKT2_2018', "sd <subdataset>" where <subdataset> is the name of a sub-dataset. See \url{https://gdal.org/programs/gdalinfo.html}}  
  \item{print}{logical. If \code{TRUE}, print the results}
  \item{open_opt}{character. Driver specific open options}   
  \item{...}{additional arguments. None implemented}   
}


\value{
character (invisibly, if \code{print=TRUE})
}


\examples{
f <- system.file("exdata/test.tif", package="terra")
gdal_info(f)

g <- gdal_info(f, c("json", "nomd", "proj4"), print=FALSE)
cat(g, "\n")
}

\keyword{spatial}

