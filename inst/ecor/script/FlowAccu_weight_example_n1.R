
library(terra)

##source("~/local/rpackages/jrc/terra/R/generics_watershed2_extension.R")


elev1 <- array(NA,c(9,9))
elev2 <- elev1

dx <- 1
dy <- 1 
for (r in 1:nrow(elev1)) {
  y <- (r-5)*dx
  for (c in 1:ncol(elev1)) {
    
    x <- (c-5)*dy
    elev1[r,c] <- 5*(x^2+y^2)
    elev2[r,c] <- 10+5*(abs(x))-0.001*y ### 5*(x^2+y^2)
  }
} 


## Elevation Raster Maps
elev1 <- rast(elev1)
elev2 <- rast(elev2)
##
weight <- elev1*0+10

t(array(elev1[],rev(dim(elev1)[1:2])))
t(array(elev2[],rev(dim(elev2)[1:2])))

plot(elev1)
plot(elev2)

## Flow Direction Raster Maps
flowdir1<- terrain(elev1,v="flowdir")
flowdir2<- terrain(elev2,v="flowdir")


t(array(flowdir1[],rev(dim(flowdir1)[1:2])))
t(array(flowdir2[],rev(dim(flowdir2)[1:2])))

plot(flowdir1)
plot(flowdir2)

## 

flow_acc1 <- flowAccu2(flowdir1)
flow_acc2 <- flowAccu2(flowdir2)

t(array(flow_acc1[],rev(dim(flow_acc1)[1:2])))
t(array(flow_acc2[],rev(dim(flow_acc2)[1:2])))

plot(flow_acc1)
plot(flow_acc2)




flow_acc1w <- flowAccu2_weight(flowdir1,weight)
flow_acc2w <- flowAccu2_weight(flowdir2,weight)

t(array(flow_acc1w[],rev(dim(flow_acc1w)[1:2])))
t(array(flow_acc2w[],rev(dim(flow_acc2w)[1:2])))

plot(flow_acc1w)
plot(flow_acc2w)


## Application to Central Africa

##library(chirps)

flowdir_f30 <- rast(system.file('ecor/data/flowdir/nigrb_dir_30s.tif',package="terra")) 
flowdir_f30[!(flowdir_f30 %in% 2^(0:7))] <- 0 

weight_f30_1 <- cellSize(flowdir_f30,unit="km")
flowacc_v_f30w_1 <- flowAccu2_weight(flowdir_f30,weight_f30_1)
flowacc_v_f30_1 <- flowAccu2(flowdir_f30)
##dates <- c("2011-01-01-15","2020-12-31")
plot(flowacc_v_f30w_1[],flowacc_v_f30_1[])
##prec <- get_chirps(flowdir_f30,dates=dates,server)
#acc_f30 <- FlowAccu_weight(flowdir_f30)
#plot(acc_f30>100)
