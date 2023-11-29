\docType{methods}
\name{watershed2}

\alias{watershed2}
\alias{watershed2,SpatRaster,integer-method}


\title{Watershed Extraction}

\description{
It extractions a watershed mask from a flow direction input raster point a pour-point cell.
}

\usage{
\S4method{watershed2}{SpatRaster,integer}(p,pp_offset,filename="",...) 
}

\arguments{
  \item{p}{Flow Direction raster map, a \code{\link{SpatRaster-class}} object, see also \code{\link{terrain}}. }
  \item{pp_offset}{pour point offset.}
  \item{filename,...}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) map containing value 1 in the extracted watershed mask and value 0 elsewhere. 
}

\details{
Lorem Ipsum. You can insert information about tha algorithm
}


\author{

Ezio Crestaz, Emanuele Cordano
}

\examples{


library(sf)




flowdir_f30 <- system.file('ecor/data/flowdir/nigrb_dir_30s.tif',package="terra") 

## path to Pourpoint shapefile:
pourpoint_f <- system.file('ecor/data/points/nigrb_pp_3s.shp',package="terra") 

###


flowdir <- rast(flowdir_f30)



pourpoint <- st_read(pourpoint_f)



offset_pp <- cellFromXY(flowdir,st_coordinates(pourpoint))
poutr <- watershed2((flowdir),as.integer(offset_pp))

\dontrun{
## Raster to Vector (polygons)
poutr[poutr==0] <- NA
poutv <- st_as_sf(as.polygons(poutr))

}




}


\references{
Lorem Ipsum
}



\keyword{spatial}
