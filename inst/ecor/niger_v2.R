
rm(list=ls())

###library(watershed2)
library(terra)
library(sf)
library(mapview)

flowdir_f <- '/stablo/local/data/hydrosheds/flowdir/af_dir_30s.bil' ## path to Flowdir Raster
pourpoint_f <- '/stablo/local/data/hydrosheds/flowdir/nigrb_pp_30s.shp' ## path to Pourpoint shapefile


flowdir <- rast(flowdir_f)

###flowdirm <- array(as.integer(flowdirm),dim(flowdirm))

pourpoint <- st_read(pourpoint_f)

##stop("here")

offset_pp <- cellFromXY(flowdir,st_coordinates(pourpoint))
###rc <- t(as.integer(dim(flowdirm)/2))


elapsed_time <- system.time({
  poutr <- watershed2((flowdir),as.integer(offset_pp))
 ## poutr <- rast(flowdir)
 ## values(poutr) <- watershed2((flowdir),as.integer(offset_pp))
 ### values(poutr) <- matrix(watershed2((flowdir),as.integer(rc[,2]-1),as.integer(rc[,1]-1)),dim(poutr),byrow=TRUE)
})