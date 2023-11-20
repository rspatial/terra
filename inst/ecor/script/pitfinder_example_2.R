rm(list=ls())

library(terra)
library(magrittr)
set.seed(34)
##elev_file <- system.file("ecor/data/small_catchment/tione_elev.tif",package="terra")

##elev <- rast(elev_file)

elev <- array(NA,c(9,9))
dx <- 1
dy <- 1 
for (r in 1:nrow(elev)) {
  y <- (r-5)*dx
  for (c in 1:ncol(elev)) {
    
    x <- (c-5)*dy
    elev[r,c] <- 10+5*(x^2+y^2) ##+1000*y
    #elev[r,c] <- y ## 10+5*(x^2+y^2)+1000*y
    }
  } 
  
elev <- cbind(elev,elev,elev,elev) 
elev <- rbind(elev,elev,elev,elev) 

elev <- rast(elev)
###elev_file <- system.file("ecor/data/small_catchment/tione_elev.tif",package="terra")

####elev <- rast(elev_file)
flowdir<- terrain(elev,v="flowdir")
t(array(flowdir[],rev(dim(flowdir)[1:2])))
###
rough <- terrain(elev,v="roughness")
tri <- terrain(elev,v="TRI")
tpi <- terrain(elev,v="TPI")

pits <- pitfinder2(flowdir)
pits <- pitfinder2(flowdir,filename="~/local/temp/pit.grd",overwrite=TRUE)
