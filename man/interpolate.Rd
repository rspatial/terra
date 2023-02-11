\name{interpolation}

\docType{methods}

\alias{interpolate}
\alias{interpolate,SpatRaster-method}

\title{Spatial interpolation}

\description{
Make a SpatRaster with interpolated values using a fitted model object of classes such as "gstat" (gstat package) or "Krige" (fields package), or any other model that has location (e.g., "x" and "y", or "longitude" and "latitude") as predictors (independent variables). If x and y are the only predictors, it is most efficient if you provide an empty (no associated data in memory or on file) SpatRaster for which you want predictions. If there are more spatial predictor variables provide these as a SpatRaster in the first argument of the function. If you do not have x and y locations as implicit predictors in your model you should use \code{\link[terra]{predict}} instead.
}

\usage{
\S4method{interpolate}{SpatRaster}(object, model, fun=predict, ..., xyNames=c("x", "y"), 
       factors=NULL, const=NULL, index = NULL, cores=1, cpkgs=NULL, 
	   na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
\item{object}{SpatRaster}
\item{model}{model object}
\item{fun}{function. Default value is "predict", but can be replaced with e.g. "predict.se" (depending on the class of \code{model}), or a custom function (see examples)}
\item{...}{additional arguments passed to \code{fun}}
\item{xyNames}{character. variable names that the model uses for the spatial coordinates. E.g., \code{c("longitude", "latitude")}}
 \item{factors}{list with levels for factor variables. The list elements should be named with names that correspond to names in \code{object} such that they can be matched. This argument may be omitted for some models from which the levels can be extracted from the \code{model} object}
\item{const}{data.frame. Can be used to add a constant for which there is no SpatRaster for model predictions. This is particularly useful if the constant is a character-like factor value}
  \item{index}{positive integer or NULL. Allows for selecting of the variable returned if the model returns multiple variables}
 \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used}
  \item{cpkgs}{character. The package(s) that need to be loaded on the nodes to be able to run the model.predict function (see examples in \code{\link{predict}})}  
 \item{na.rm}{logical. If \code{TRUE}, cells with \code{NA} values in the predictors are removed from the computation. This option prevents errors with models that cannot handle \code{NA} values. In most other cases this will not affect the output. An exception is when predicting with a model that returns predicted values even if some (or all!) variables are \code{NA} }

\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{predict}}, \code{link{interpIDW}}, \code{link{interpNear}}}


\examples{

r <- rast(system.file("ex/elev.tif", package="terra"))
ra <- aggregate(r, 10)
xy <- data.frame(xyFromCell(ra, 1:ncell(ra)))
v <- values(ra)
i <- !is.na(v)
xy <- xy[i,]
v <- v[i]

\dontrun{
library(fields) 
tps <- Tps(xy, v)
p <- rast(r)

# use model to predict values at all locations
p <- interpolate(p, tps)
p <- mask(p, r)
plot(p)

### change "fun" from predict to fields::predictSE to get the TPS standard error
## need to use "rast(p)" to remove the values
se <- interpolate(rast(p), tps, fun=predictSE)
se <- mask(se, r)
plot(se)

### another predictor variable, "e"
e <- (init(r, "x") * init(r, "y")) / 100000000
names(e) <- "e"

z <- as.matrix(extract(e, xy)[,-1])

## add as another independent variable
xyz <- cbind(xy, z)
tps2 <- Tps(xyz, v)
p2 <- interpolate(e, tps2, xyOnly=FALSE)

## as a linear coveriate
tps3 <- Tps(xy, v, Z=z)

## Z is a separate argument in Krig.predict, so we need a new function
## Internally (in interpolate) a matrix is formed of x, y, and elev (Z)

pfun <- function(model, x, ...) {
   predict(model, x[,1:2], Z=x[,3], ...)
}
p3 <- interpolate(e, tps3, fun=pfun)


#### gstat examples
library(gstat)
library(sp)
data(meuse)

### inverse distance weighted (IDW)
r <- rast(system.file("ex/meuse.tif", package="terra"))
mg <- gstat(id = "zinc", formula = zinc~1, locations = ~x+y, data=meuse, 
            nmax=7, set=list(idp = .5))
z <- interpolate(r, mg, debug.level=0, index=1)
z <- mask(z, r)

## with a model built with an `sf` object you need to provide custom function

library(sf)
sfmeuse <- st_as_sf(meuse, coords = c("x", "y"), crs=crs(r))
mgsf <- gstat(id = "zinc", formula = zinc~1, data=sfmeuse,  nmax=7, set=list(idp = .5))

interpolate_gstat <- function(model, x, crs, ...) {
	v <- st_as_sf(x, coords=c("x", "y"), crs=crs)
	p <- predict(model, v, ...)
	as.data.frame(p)[,1:2]
}

zsf <- interpolate(r, mgsf, debug.level=0, fun=interpolate_gstat, crs=crs(r), index=1)
zsf <- mask(zsf, r)

### kriging

### ordinary kriging
v <- variogram(log(zinc)~1, ~x+y, data=meuse)
mv <- fit.variogram(v, vgm(1, "Sph", 300, 1))
gOK <- gstat(NULL, "log.zinc", log(zinc)~1, meuse, locations=~x+y, model=mv)
OK <- interpolate(r, gOK, debug.level=0)

## universal kriging
vu <- variogram(log(zinc)~elev, ~x+y, data=meuse)
mu <- fit.variogram(vu, vgm(1, "Sph", 300, 1))
gUK <- gstat(NULL, "log.zinc", log(zinc)~elev, meuse, locations=~x+y, model=mu)
names(r) <- "elev"
UK <- interpolate(r, gUK, debug.level=0)

## co-kriging
gCoK <- gstat(NULL, 'log.zinc', log(zinc)~1, meuse, locations=~x+y)
gCoK <- gstat(gCoK, 'elev', elev~1, meuse, locations=~x+y)
gCoK <- gstat(gCoK, 'cadmium', cadmium~1, meuse, locations=~x+y)
gCoK <- gstat(gCoK, 'copper', copper~1, meuse, locations=~x+y)
coV <- variogram(gCoK)
plot(coV, type='b', main='Co-variogram')
coV.fit <- fit.lmc(coV, gCoK, vgm(model='Sph', range=1000))
coV.fit
plot(coV, coV.fit, main='Fitted Co-variogram')
coK <- interpolate(r, coV.fit, debug.level=0)
plot(coK)
}
}

\keyword{methods}
\keyword{spatial}
