rm(list=ls())
library(terra)

library(calculus)

color_palette <- terrain.colors(50)

# 32	64	128
# 16	x	  1
# 8	  4	  2
flowdir2degN <- function(x,clockwise=TRUE,angle_from_x_axis=90){
  out <- -log(x)/log(2)*45
  out <-  out-angle_from_x_axis
  if (clockwise) out <- -out
  return(out)
} ## see help(terrain)
  


### ---------------------------------------------------------
### 1. GRID DEFINITION
### ---------------------------------------------------------
dx <- dy <- 1
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

model=c("vshape","paraboloid","planar")[1]

if (model=="vshape") {
  elev_f1 <-  "2.0*x*tanh(100*x) + 0.5*y"  
  elev_f2 <-   "2.0*x*tanh(100*x) - 0.5*y"
} else if ((model=="paraboloid")) {
  elev_f1 <- " 10*(x^2+y^2)^0.8"
  elev_f2 <- " -10*(x^2+y^2)^0.8"
} else if ((model=="planar")) {
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

names(r1) <- "elev"
names(r2) <- "elev"

### ---------------------------------------------------------
### 6. FLOW DIRECTIONS (terra + LAD/LTD)
### ---------------------------------------------------------
deviation_types=c("lad","ltd")
lambdas <- c(0,0.5,1.0)

# Elevation Model A
r1$flow_d8 <- terrain(r1$elev, "flowdir")
r1$aspect <-  terrain(r1$elev, "aspect")
# Elevation Model B
r2$flow_d8 <- terrain(r2$elev, "flowdir")
r2$aspect <-  terrain(r2$elev, "aspect")
###
dev1 <- list()
cumdev1 <- list()
cumavgdev <- list()
dev2 <- list()
cumdev2 <- list()
####
for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    it <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    r1[[it]] <- flowDir(r1$elev, lambda=lambda, deviation_type=deviation_type)
    r2[[it]] <- flowDir(r2$elev, lambda=lambda, deviation_type=deviation_type)
    ####
    dev1[[it]] <- flowdir2degN(r1[[it]])-r1$aspect
    dd <- dev1[[it]]
    dd[is.na(dd)] <- 0
    cumdev1[[it]] <- flowAccumulation(r1[[it]],weight=dd)
    cumdev1[[it]][is.na(dev1[[it]])] <- NA 
    cumdev1[[it]][NIDP(r1[[it]])!=1] <- NA 
    faf <- flowAccumulation(r1[[it]])
    cumavgdev[[it]] <- cumdev1[[it]]/faf
    
    
    
  }
}

### ---------------------------------------------------------
### 7. VISUALIZATION EXAMPLES
### ---------------------------------------------------------
nrr <- 1+length(lambdas)*length(deviation_types)
par(mfrow=c(nrr,2))

# Model A
plot(r1$elev, main="Elevation Model A", col=color_palette)
terra::as.arrows(r1$flow_d8, unit="flowdir", col="black")
as.arrows(r1$aspect, unit="degrees", col="blue",clockwise=TRUE,angle_from_x_axis=90)
# Model B
plot(r2$elev, main="Elevation Model B", col=color_palette)
terra::as.arrows(r2$flow_d8, unit="flowdir", col="black")
terra::as.arrows(r2$aspect, unit="degrees", col="blue",clockwise=TRUE,angle_from_x_axis=90)
for (lambda in lambdas) {
  for (deviation_type in deviation_types) {
    ita <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    plot(r1$elev, main=sprintf("Model A - %s λ=%1.2f",toupper(deviation_type),lambda),
         col=color_palette)
    terra::as.arrows(r1[[ita]], unit="flowdir", col="black")
    terra::as.arrows(r1$aspect, unit="degrees", col="blue",clockwise=TRUE,angle_from_x_axis=90)
    itb <- sprintf("flow_%s%03d", deviation_type, round(lambda*10))
    plot(r2$elev, main=sprintf("Model B - %s λ=%1.2f",toupper(deviation_type),lambda),
         col=color_palette)
    terra::as.arrows(r2[[itb]], unit="flowdir", col="black")
    terra::as.arrows(r2$aspect, unit="degrees", col="blue",clockwise=TRUE,angle_from_x_axis=90)
  }
}

par(mfrow=c(1,1))

message("Analysis complete")
# plot(cumavgdev$flow_ltd010)
# > plot(cumavgdev$flow_lad010)
# > plot(cumavgdev$flow_lad000)
# > plot(cumavgdev$flow_ltd000)
# > plot(cumavgdev$flow_lad005)
# > plot(cumavgdev$flow_ltd005)
# > source("~/local/rpackages/jrc/terra/inst/ltd_flow_tests/flowdir_example_vv2_with_aspect.R")
# Analysis complete
# > plot(cumavgdev$flow_ltd005)
# > plot(cumavgdev$flow_ltd000)
# > plot(cumavgdev$flow_ltd000<5)
# > plot(abs(cumavgdev$flow_ltd000)<5)
# > plot(abs(cumavgdev$flow_ltd000)<10)
# > plot(abs(cumavgdev$flow_ltd000)<20)
# > plot(abs(cumavgdev$flow_ltd000)<15)
# > plot(abs(cumavgdev$flow_ltd000)<10)
# > plot(abs(cumavgdev$flow_ltd000)<12)
# > plot(abs(cumavgdev$flow_ltd000)<11)
# > plot(abs(cumavgdev$flow_ltd000)<13)
# > plot(abs(cumavgdev$flow_ltd000)<12.5)
# > plot(abs(cumavgdev$flow_ltd050)<12.5)
# Error in h(simpleError(msg, call)) : 
#   error in evaluating the argument 'x' in selecting a method for function 'plot': non-numeric argument to mathematical function
# > plot(abs(cumavgdev$flow_ltd010)<12.5)
# > plot(abs(cumavgdev$flow_ltd000)<12.5)
# > plot(abs(cumavgdev$flow_ltd010)<10)
# > plot(abs(cumavgdev$flow_ltd010)<5)
# 
# 
# 
