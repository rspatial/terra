rm(list=ls())

library(calculus)
library(terra)
library(stringr)

### ---------------------------------------------------------
### 1. GRID DEFINITION
### ---------------------------------------------------------
options(terra.pal=terrain.colors(10))
dx <- dy <- 1 ## 10
xmin <- -10
xmax <- 10
ymin <- -10
ymax <- 10

x <- seq(xmin, xmax, by = dx)
y <- seq(ymin, ymax, by = dy)

vars <- expand.grid(x = x, y = y)

### ---------------------------------------------------------
### 2. DEFINE TWO ANALYTICAL ELEVATION FUNCTIONS
### ---------------------------------------------------------


model=c("vshape","paraboloid","planar")[1] ## set vshape (if 1) or planar (if 3) for instance 

if (model=="vshape") {
  # Model: V-shape valley with longitudinal slope
  elev_f1 <-  "2.0*x*tanh(100*x) + 0.5*y"  
  elev_f2 <-   "2.0*x*tanh(100*x) - 0.5*y"
} else if ((model=="paraboloid")) {
 # Model: Parabolide longitudinal slope
 elev_f1 <- " 10*(x^2+y^2)^0.8"
 elev_f2 <- " -10*(x^2+y^2)^0.8"
} else if ((model=="planar")) {
  # Model: Parabolic longitudinal slope
  elev_f1 <- " 2.0*x+0.225*y"
  elev_f2 <- "2.0*x-0.225*y"
}else {
  
  elev_f1 <- NULL
  elev_f2 <- NULL
  
}
### ---------------------------------------------------------
### 3. EVALUATE BOTH ELEVATION MODELS
### ---------------------------------------------------------

vars$elev1 <- evaluate(elev_f1, var = vars[,c("x","y")],params=list(ymin=ymin,ymax=ymax))
vars$elev2 <- evaluate(elev_f2, var = vars[,c("x","y")],params=list(ymin=ymin,ymax=ymax))



### ---------------------------------------------------------
### 5. BUILD RASTERS
### ---------------------------------------------------------

r1 <- rast(vars[,c("x","y","elev1")], type="xyz")
r2 <- rast(vars[,c("x","y","elev2")], type="xyz")
###
names(r1) <- "elev"
names(r2) <- "elev"
####


### ---------------------------------------------------------
### 6. FLOW DIRECTIONS (terra + LAD/LTD)
### ---------------------------------------------------------
deviation_types=c("lad","ltd") ## set LTD for instance
lambdas <- c(0,0.5,1)

# Elevation Model A
r1$flow_d8      <- terrain(r1$elev, "flowdir")
# Elevation Model B 
r2$flow_d8      <- terrain(r1$elev, "flowdir")
for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    it <- "flow_%s%03d"|> sprintf(deviation_type,round(lambda*10))
    # Elavation Model A
    r1[[it]] <- flowDir(r1$elev, lambda=lambda,   deviation_type=deviation_type)
    # Elevation Model B 
    r2[[it]] <- flowDir(r2$elev, lambda=lambda,   deviation_type=deviation_type)
  }
  
  
}

### ---------------------------------------------------------
### 7. VISUALIZATION EXAMPLES
### ---------------------------------------------------------
nrr <- 1+length(lambdas)*length(deviation_types)
par(mfrow=c(nrr,2))

# Model A
plot(r1$elev, main="Elevation Model A")
as.arrows(r1$flow_d8, unit="flowdir", col="black")
# Model B
plot(r2$elev, main=sprintf("Elevation Model B ",toupper(deviation_type)))
as.arrows(r2$flow_d8, unit="flowdir", col="black")
for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    ita <- "flow_%s%03d"|> sprintf(deviation_type,round(lambda*10))
    # Elavation Model A
    plot(r1$elev, main=sprintf("Model A - %s λ=%1.2f",toupper(deviation_type),lambda))
    terra::as.arrows(r1[[ita]], unit="flowdir", col="black")
    # Elevation Model B 
    itb <- "flow_%s%03d"|> sprintf(deviation_type,round(lambda*10))
    plot(r2$elev, main=sprintf("Model B - %s λ=%1.2f",toupper(deviation_type),lambda))
    terra::as.arrows(r2[[itb]], unit="flowdir", col="black")
  }
  
}



par(mfrow=c(1,1))

message("Analysis complete")

