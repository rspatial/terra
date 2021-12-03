#
# Testing "watershed2" Niger River Basin 
#
rm(list=ls())

library(terra)
library(sf)
library(mapview)

flowdir_f <- '/stablo/local/data/hydrosheds/flowdir/af_dir_30s.bil' ## path to Flowdir Raster
poutr_f <- '/stablo/local/data/hydrosheds/flowdir/niger_basin.tif'
pourpoint_f <- '/stablo/local/data/hydrosheds/flowdir/nigrb_pp_30s.shp' ## path to Pourpoint shapefile

flowdir <- rast(flowdir_f)
pourpoint <- st_read(pourpoint_f)


offset_pp <- cellFromXY(flowdir,st_coordinates(pourpoint))

## S4 Method 'watershed2' 

## Uncomment the following line to see the "watershed2" documentation
# help(watershed2)


elapsed_time <- system.time({
  poutr <- watershed2((flowdir),as.integer(offset_pp),filename=poutr_f,overwrite=TRUE)

})

