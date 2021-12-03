
rm(list=ls())

###library(watershed2)
library(terra)
library(sf)
library(mapview)

##flowdir_f <- '/stablo/local/data/hydrosheds/flowdir/af_dir_30s.bil' ## path to Flowdir Raster
##pourpoint_f <- '/stablo/local/data/hydrosheds/flowdir/PourPoint_wgs84_Niger.shp' ## path to Pourpoint shapefile
set.seed(22)

elev <- rast('/home/ecor/local/rpackages/terra_rspatial/inst/ecor/input/dem.asc')
telev <- 9999
elev[is.na(elev)] <- telev

flowdir <- terrain(elev,"flowdir")
flowdir[elev==telev] <- 247


##flowdir <- writeRaster(flowdir,filename=paste0(tempfile(),".grd"),overwrite=TRUE)
##stop("CONTROLLA GLI NANs!!!")
###pourpoint <- st_coordinetes(st_read(pourpoint_f))
pourpoint <- as.matrix(t(c(x=69,y=15)))
pp_offset <- cellFromXY(flowdir,pourpoint)

pourpoint <- as.matrix(t(c(x=60,y=30)))
##pp_offset <- cellFromXY(flowdir,pourpoint)
###flowdir[pp_offset] <- -10

###rc <- t(as.integer(dim(flowdirm)/2))


elapsed_time <- system.time({
  poutr <- watershed2(flowdir,as.integer(pp_offset)) ###,filename=paste0(tempfile(),".grd"),overwrite=TRUE)
 ####   matrix(watershed2((flowdir),as.integer(rc[,2]-1),as.integer(rc[,1]-1)),dim(flowdir),byrow=TRUE)
})