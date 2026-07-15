library(terra)
options(terra.pal=terrain.colors(10))
rr <- list()
rr$elev <- rast(system.file('ex/elev_vinschgau.tif', package="terra"))

### ---------------------------------------------------------
### 6. FLOW DIRECTIONS (terra + LAD/LTD)
### ---------------------------------------------------------
deviation_types=c("lad","ltd") ## set LTD for instance
lambdas <- c(0,0.25,0.5,0.75,1)

# Elevation Model A
rr$flow_d8      <- terrain(rr$elev, "flowdir")

for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    it <- "flow_%s%03d"|> sprintf(deviation_type,round(lambda*10))
    rr[[it]] <- flowdirD8ltd(rr$elev, lambda=lambda,   deviation_type=deviation_type)
   
  }
  
  
}
nrr <- 1+length(lambdas)*length(deviation_types)
par(mfrow=c(nrr,2))

plot(rr$elev, main=sprintf("Elevation Model B ",toupper(deviation_type)))
arrows_on_rast(rr$flow_d8, unit="flowdir", col="black",arrow_length = 0.001)
for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    it <- "flow_%s%03d"|> sprintf(deviation_type,round(lambda*10))
  
    plot(rr$elev, main=sprintf("Model A - %s λ=%1.2f",toupper(deviation_type),lambda))
    arrows_on_rast(rr[[it]], unit="flowdir", col="black",arrow_length = 0.001)
    
  }
  
}


par(mfrow=c(1,1))

message("Analysis complete")
