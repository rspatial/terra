\name{interpolate}

\docType{methods}

\alias{interpolate}
\alias{interpolate,SpatRaster-method}

\title{Interpolate}

\description{
Make a RasterLayer with interpolated values using a fitted model object of classes such as "gstat" (gstat package) or "Krige" (fields package). That is, these are models that have location ("x" and "y", or "longitude" and "latitude") as independent variables. If x and y are the only independent variables provide an empty (no associated data in memory or on file) SpatRaster for which you want predictions. If there are more spatial predictor variables provide these as a Raster* object in the first argument of the function. If you do not have x and y locations as implicit predictors in your model you should use \code{\link[terra]{predict}} instead.
}

\usage{
\S4method{interpolate}{SpatRaster}(object, model, fun=predict, ..., 
       xyNames=c("x", "y"), factors=NULL, const=NULL, index = NULL, 
	   na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
\item{object}{SpatRaster}
\item{model}{model object}
\item{fun}{function. Default value is "predict", but can be replaced with e.g. "predict.se" (depending on the class of the model object)}
\item{...}{additional arguments passed to \code{fun}}
\item{xyNames}{character. variable names that the model uses for the spatial coordinates. E.g., \code{c("longitude", "latitude")}}
 \item{factors}{list with levels for factor variables. The list elements should be named with names that correspond to names in \code{object} such that they can be matched. This argument may be omitted for many models as the predict function will extract the levels from the \code{model} object}
\item{const}{data.frame. Can be used to add a constant for which there is no Raster object for model predictions. This is particularly useful if the constant is a character-like factor value}
  \item{index}{positive integer or NULL. Allows for selecting of the variable returned if the model returns multiple variables}
 \item{na.rm}{logical. If \code{TRUE}, cells with \code{NA} values in the predictors are removed from the computation. This option prevents errors with models that cannot handle \code{NA} values. In most other cases this will not affect the output. An exception is when predicting with a model that returns predicted values even if some (or all!) variables are \code{NA} }

\item{filename}{character. Output filename. Optional}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{predict}}}



\examples{
\donttest{
## Thin plate spline interpolation with x and y only
# some example data
r <-rast(system.file("exdata/test.tif", package="terra"))
ra <- aggregate(r, 10)
xy <- data.frame(xyFromCell(ra, 1:ncell(ra)))
v <- values(ra)

#### Thin plate spline model
library(fields) 
tps <- Tps(xy, v)
x <- rast(r)

# use model to predict values at all locations
p <- interpolate(x, tps)
p <- mask(p, r)
plot(p)

## change the fun from predict to fields::predictSE to get the TPS standard error
se <- interpolate(x, tps, fun=predictSE)
se <- mask(se, r)
plot(se)

## another variable; let"s call it elevation
elevation <- (init(r, "x") * init(r, "y")) / 100000000
names(elevation) <- "elev"
elevation <- mask(elevation, r)

z <- extract(elevation, vect(xy), fun=function(x)x[1])

# add as another independent variable
vv <- na.omit(cbind(xy, z, v))
tps2 <- Tps(vv[,1:3], vv[,4])
p2 <- interpolate(elevation, tps2)
plot(p2) 

# as a linear coveriate
tps3 <- Tps(vv[,1:2], vv[,4], Z=vv[,3])

# Z is a separate argument in Krig.predict, so we need a new function
# Internally (in interpolate) a matrix is formed of x, y, and elev (Z)

pfun <- function(model, x, ...) {
   predict(model, x[,1:2], Z=x[,3], ...)
}
p3 <- interpolate(elevation, tps3, fun=pfun)

#### gstat examples
library(sp)
library(gstat)
data(meuse)

## inverse distance weighted (IDW)
r <-rast(system.file("exdata/test.tif", package="terra"))
mg <- gstat(id = "zinc", formula = zinc~1, locations = ~x+y, data=meuse, 
            nmax=7, set=list(idp = .5))
f <- function(model, data, ...) predict(model, data, ...)[,3,drop=FALSE]
z <- interpolate(r, mg, fun=f, debug.level=0)

## kriging
coordinates(meuse) <- ~x+y
crs(meuse) <- crs(r)

## ordinary kriging
v <- variogram(log(zinc)~1, meuse)
m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
gOK <- gstat(NULL, "log.zinc", log(zinc)~1, meuse, model=m)

fg <- function(model, d, crs, ...) {
	sp <- SpatialPointsDataFrame(d[,1:2,drop=FALSE], data.frame(d), proj4string=CRS(crs))
	data.frame(predict(model, sp, ...))[,3:4]
}
OK1 <- interpolate(r, gOK, fun=fg, debug.level=0, crs=crs(r), na.rm=TRUE)
OK2 <- interpolate(r, gOK, fun=fg, debug.level=0, crs=crs(r), na.rm=FALSE)
OK3 <- interpolate(x, gOK, fun=fg, debug.level=0, crs=crs(x))
plot(c(OK1[[1]], OK2[[1]], OK3[[1]]))


# examples below provided by Maurizio Marchi
## universial kriging
vu <- variogram(log(zinc)~elev, meuse)
mu <- fit.variogram(vu, vgm(1, "Sph", 300, 1))
gUK <- gstat(NULL, "log.zinc", log(zinc)~elev, meuse, model=mu)
names(r) <- "elev"
UK <- interpolate(r, gUK, fun=fg, debug.level=0, crs=crs(r))

## co-kriging
gCoK <- gstat(NULL, "log.zinc", log(zinc)~1, meuse)
gCoK <- gstat(gCoK, "elev", elev~1, meuse)
gCoK <- gstat(gCoK, "cadmium", cadmium~1, meuse)
gCoK <- gstat(gCoK, "copper", copper~1, meuse)
coV <- variogram(gCoK)
plot(coV, type="b", main="Co-variogram")
coV.fit <- fit.lmc(coV, gCoK, vgm(model="Sph", range=1000))
coV.fit
plot(coV, coV.fit, main="Fitted Co-variogram")
coK <- interpolate(r, coV.fit, fun=fg, debug.level=0, crs=crs(r))
plot(coK)
}
}

\keyword{methods}
\keyword{spatial}
