
rm(list=ls())

###library(watershed2)
library(terra)
library(sf)
library(mapview)


flowdir_f3 <- '/home/ecor/local/data/spatial/hydrosheds/flowdir/niger/nigrb_dir_3s.bil' ###af_dir_30s.bil' ## path to Flowdir Raster
flowdir_f30 <- '/home/ecor/local/rpackages/jrc/terra/inst/ecor/data/flowdir/nigrb_dir_30s.tif'   ##data/spatial/hydrosheds/flowdir/niger/nigrb_dir_30s.tif'



flowdir_f <- '/home/ecor/local/data/spatial/hydrosheds/flowdir/af_dir_30s.bil' ##poutrniger/nigrb_dir_3s.bil' ###af_dir_30s.bil' ## path to Flowdir Raster
pourpoint_f <- '/home/ecor/local/data/spatial/hydrosheds/flowdir/niger/nigrb_pp_3s.shp' ## path to Pourpoint shapefile
pourpoint_30s_f <- '/home/ecor/local/data/spatial/hydrosheds/flowdir/PourPoint_wgs84_Niger.shp' ##   nigrb_pp_30s.shp' ## path to Pourpoint shapefile

###


###


flowdir <- rast(flowdir_f)

flowdir_v2 <- crop(flowdir,rast(flowdir_f3),filename=flowdir_f30,overwrite=TRUE)
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