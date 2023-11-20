

library(terra)


flowdir_f30 <- system.file('ecor/data/flowdir/nigrb_dir_30s.tif',package="terra") 


flowdir <- rast(flowdir_f30)


poutr <- NIDP2((flowdir))

#########

elev <- array(NA,c(9,9))
dx <- 1
dy <- 1 
for (r in 1:nrow(elev)) {
  y <- (r-5)*dx
  for (c in 1:ncol(elev)) {
    
    x <- (c-5)*dy
    elev[r,c] <- 10+5*(abs(x))-0.001*y ### 5*(x^2+y^2)
  }
} 


#elev <- array(c(1000,1000,11,11,11,10.5,10.5,10,0),c(3,3))


elev <- rast(elev)
flowdir<- terrain(elev,v="flowdir")
t(array(flowdir[],rev(dim(flowdir)[1:2])))
pit <- pitfinder2((flowdir))
nidp <- NIDP2((flowdir))
t(array(nidp[],rev(dim(flowdir)[1:2])))
