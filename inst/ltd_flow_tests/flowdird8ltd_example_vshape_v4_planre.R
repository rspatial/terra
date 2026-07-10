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


model=c("vshape","paraboloid","planar")[1] 

if (model=="vshape") {
  # Model: V-shape valley with longitudinal slope
  elev_f1 <-  "2.0*x*tanh(100*x) + 0.5*y"  ##"x^2+y^2" ##
  elev_f2 <-   "2.0*x*tanh(100*x) - 0.5*y"
} else if ((model=="paraboloid")) {
 # Model: Parabolide longitudinal slope
 elev_f1 <- " 10*(x^2+y^2)^0.8"
 elev_f2 <- " -10*(x^2+y^2)^0.8"
} else if ((model=="planar")) {
  # Model: Parabolide longitudinal slope
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

# Remove borders (optional)
#border_idx <- which(vars$x %in% range(vars$x) | vars$y %in% range(vars$y))
#vars$elev1[border_idx] <- NA
#vars$elev2[border_idx] <- NA

### ---------------------------------------------------------
### 4. ANALYTICAL GRADIENTS (Model A & B)
### ---------------------------------------------------------

# grad1 <- matrix(gradient(elev_f1, vars[,c("x","y")], accuracy = 8), ncol = 2)
# grad2 <- matrix(gradient(elev_f2, vars[,c("x","y")], accuracy = 8), ncol = 2)
# 
# vars$flow_angle1 <- atan2(-grad1[,2], -grad1[,1])
# vars$flow_angle2 <- atan2(-grad2[,2], -grad2[,1])

### ---------------------------------------------------------
### 5. BUILD RASTERS
### ---------------------------------------------------------

r1 <- rast(vars[,c("x","y","elev1")], type="xyz")
r2 <- rast(vars[,c("x","y","elev2")], type="xyz")

### ---------------------------------------------------------
### 6. FLOW DIRECTIONS (terra + LAD/LTD)
### ---------------------------------------------------------
deviation_type="lad"
# Model A

r1$flow_d8      <- terrain(r1$elev1, "flowdir")
r1$flow_lad0    <- flowdirD8ltd(r1$elev1, lambda=0,   deviation_type=deviation_type)

r1$flow_lad05   <- flowdirD8ltd(r1$elev1, lambda=0.75, deviation_type=deviation_type)
r1$flow_lad1    <- flowdirD8ltd(r1$elev1, lambda=1, deviation_type=deviation_type)
# Model B
r2$flow_d8      <- terrain(r2$elev2, "flowdir")
r2$flow_lad0    <- flowdirD8ltd(r2$elev2, lambda=0,   deviation_type=deviation_type)
r2$flow_lad05   <- flowdirD8ltd(r2$elev2, lambda=0.75, deviation_type=deviation_type)
r2$flow_lad1    <- flowdirD8ltd(r2$elev2, lambda=1,   deviation_type=deviation_type)

### ---------------------------------------------------------
### 7. VISUALIZATION EXAMPLES
### ---------------------------------------------------------

par(mfcol=c(4,2))

# Model A
plot(r1$elev1, main="Elevation Model A")
arrows_on_rast(r1$flow_d8, unit="flowdir", col="black")

plot(r1$elev1, main=sprintf("Model A - %s λ=0",toupper(deviation_type)))
arrows_on_rast(r1$flow_lad0, unit="flowdir", col="black")


plot(r1$elev1, main=sprintf("Model A - %s λ=0.5",toupper(deviation_type)))
arrows_on_rast(r1$flow_lad05, unit="flowdir", col="black")


plot(r1$elev1, main=sprintf("Model A - %s λ=1",toupper(deviation_type)))
arrows_on_rast(r1$flow_lad1, unit="flowdir", col="black")

# Model B
plot(r2$elev2, main=sprintf("Elevation Model B ",toupper(deviation_type)))
arrows_on_rast(r2$flow_d8, unit="flowdir", col="black")

plot(r2$elev2, main=sprintf("Model B - %s λ=0",toupper(deviation_type)))
arrows_on_rast(r2$flow_lad0, unit="flowdir", col="black")

plot(r2$elev2, main=sprintf("Model B - %s λ=0.5",toupper(deviation_type)))
arrows_on_rast(r2$flow_lad05, unit="flowdir", col="black")

plot(r2$elev2, main=sprintf("Model B - %s λ=1",toupper(deviation_type)))
arrows_on_rast(r2$flow_lad1, unit="flowdir", col="black")


par(mfrow=c(1,1))

##stop("Analysis complete")

