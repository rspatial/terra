
rm(list=ls())

###library(watershed2)
library(terra)
library(sf)
library(mapview)

##flowdir_f <- '/stablo/local/data/hydrosheds/flowdir/af_dir_30s.bil' ## path to Flowdir Raster
##pourpoint_f <- '/stablo/local/data/hydrosheds/flowdir/PourPoint_wgs84_Niger.shp' ## path to Pourpoint shapefile

elev <- rast('/home/ecor/local/rpackages/terra_rspatial/inst/ecor/input/dem.asc')
elev[is.na(elev)] <- -10

flowdir <- terrain(elev,"flowdir")

##stop("CONTROLLA GLI NANs!!!")
###pourpoint <- st_coordinetes(st_read(pourpoint_f))
pourpoint <- as.matrix(t(c(x=70,y=15)))
pp_offset <- cellFromXY(flowdir,pourpoint)

pourpoint <- as.matrix(t(c(x=60,y=30)))
pp_offset <- cellFromXY(flowdir,pourpoint)
##flowdir[pp_offset] <- -10

###rc <- t(as.integer(dim(flowdirm)/2))


elapsed_time <- system.time({
  poutr <- rast(flowdir)
  values(poutr) <- watershed2(flowdir,as.integer(pp_offset))
 ####   matrix(watershed2((flowdir),as.integer(rc[,2]-1),as.integer(rc[,1]-1)),dim(flowdir),byrow=TRUE)
})