rm(list=ls())

library(terra)
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
## NDPI  Raster Maps
nidp1 <- NIDP2((flowdir1))
nidp2 <- NIDP2((flowdir2))

t(array(nidp1[],rev(dim(nidp1)[1:2])))
t(array(nidp2[],rev(dim(nidp2)[1:2])))



## 
flow_acc1 <- flowAccu2((flowdir1))
flow_acc2 <- flowAccu2((flowdir2))

t(array(flow_acc1[],rev(dim(flow_acc1)[1:2])))
t(array(flow_acc2[],rev(dim(flow_acc2)[1:2])))

plot(flow_acc1)
plot(flow_acc2)

##########

flowdir_f30 <- rast(system.file('ecor/data/flowdir/nigrb_dir_30s.tif',package="terra")) 
acc_f30 <- flowAccu2(flowdir_f30)
plot(acc_f30>100)

