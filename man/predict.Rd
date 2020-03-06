\name{predict}

\docType{methods}

\alias{predict}
\alias{predict,SpatRaster-method}

\title{Spatial model predictions}

\description{
Make a SpatRaster object with predictions from a fitted model object (for example, obtained with \code{glm} or \code{randomForest}). The first argument is a SpatRaster object with the predictor variables. The \code{\link{names}} in the Raster object should exactly match those expected by the model. Any regression like model for which a predict method has been implemented (or can be implemented) can be used. 

This approach of using model predictions is commonly used in remote sensing (for the classification of satellite images) and in ecology, for species distribution modeling.
}

\usage{
\S4method{predict}{SpatRaster}(object, model, fun=predict, ..., factors=NULL, 
    const=NULL, na.rm=FALSE, index=NULL, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{object}{SpatRaster}
  \item{model}{fitted model of any class that has a "predict" method (or for which you can supply a similar method as \code{fun} argument. E.g. glm, gam, or randomForest}
  \item{fun}{function. Default value is \code{predict}, but can be replaced with e.g. predict.se (depending on the type of model), or your own custom function}
  \item{...}{additional arguments for \code{fun}}    

  \item{const}{data.frame. Can be used to add a constant value as a predictor variable so that you do not need to make a SpatRaster layer for it}
 \item{factors}{list with levels for factor variables. The list elements should be named with names that correspond to names in \code{object} such that they can be matched. This argument may be omitted for standard models such as "glm" as the predict function will extract the levels from the \code{model} object, but it is necessary in some other cases (e.g. cforest models from the party package)}
  \item{na.rm}{logical. If \code{TRUE}, cells with \code{NA} values in the predictors are removed from the computation. This option prevents errors with models that cannot handle \code{NA} values. In most other cases this will not affect the output. An exception is when predicting with a model that returns predicted values even if some (or all!) variables are \code{NA} }
  \item{index}{integer. To select subset of output variables}
 
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
 }


\value{
SpatRaster
}


\examples{
logo <- rast(system.file("exdata/logo.tif", package="terra"))   
names(logo) <- c("red", "green", "blue")
p <- matrix(c(48, 48, 48, 53, 50, 46, 54, 70, 84, 85, 74, 84, 95, 85, 
   66, 42, 26, 4, 19, 17, 7, 14, 26, 29, 39, 45, 51, 56, 46, 38, 31, 
   22, 34, 60, 70, 73, 63, 46, 43, 28), ncol=2)

a <- matrix(c(22, 33, 64, 85, 92, 94, 59, 27, 30, 64, 60, 33, 31, 9,
   99, 67, 15, 5, 4, 30, 8, 37, 42, 27, 19, 69, 60, 73, 3, 5, 21,
   37, 52, 70, 74, 9, 13, 4, 17, 47), ncol=2)

xy <- rbind(cbind(1, p), cbind(0, a))

# extract predictor values for points
e <- extract(logo, xy[,2:3])

# combine with response 
v <- data.frame(cbind(pa=xy[,1], e))

#build a model, here with glm 
model <- glm(formula=pa~., data=v)

#predict to a raster
r1 <- predict(logo, model)

plot(r1)
points(p, bg='blue', pch=21)
points(a, bg='red', pch=21)

# logistic regression
model <- glm(formula=pa~., data=v, family="binomial")
r1log <- predict(logo, model, type="response")

# use a modified function to get the probability and standard error
# from the glm model. The values returned by "predict" are in a list,
# and this list needs to be transformed to a matrix

predfun <- function(model, data) {
  v <- predict(model, data, se.fit=TRUE)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}

r2 <- predict(logo, model, fun=predfun)

# principal components of a SpatRaster
# here using sampling to simulate an object too large
# to feed all its values to prcomp

sr <- values(spatSample(logo, 100, as.raster=TRUE))
pca <- prcomp(sr)

x <- predict(logo, pca)
plot(x)

\donttest{
library(pls)
model <- plsr(formula=pa~., data=v)
# this returns an array:
predict(model, v[1:5,])
# write a function to turn that into a matrix
pfun <- function(x, data) {
   y <- predict(x, data)
   d <- dim(y)
   dim(y) <- c(prod(d[1:2]), d[3])
   y
}

pp <- predict(logo, model, fun=pfun)

# Random Forest

library(randomForest)
rfmod <- randomForest(pa ~., data=v)

## note the additional argument "type='response'" that is 
## passed to predict.randomForest
r3 <- predict(logo, rfmod, type='response')

## get class membership probabilities
vv <- v
vv$pa <- as.factor(vv$pa)
rfmod2 <- randomForest(pa ~., data=vv)
r4 <- predict(logo, rfmod2, type='prob')
rr4 <- as(r4, "Raster")
raster::spplot(rr4)


# cforest (other Random Forest implementation) example with factors argument

library(party)
m <- cforest(pa~., control=cforest_unbiased(mtry=3), data=v)
# the second argument in party:::predict.RandomForest
# is "OOB", and not "newdata" or similar. We need to write a wrapper
# predict function to deal with this 	
predfun <- function(m, d, ...) predict(m, newdata=d, ...)

pc <- predict(logo, m, OOB=TRUE, fun=predfun)

# knn, using "app" instead of "predict"
library(class)
cl <- factor(c(rep(1, nrow(p)), rep(0, nrow(a))))
train <- extract(logo, rbind(p, a))
k <- app(logo, function(x) as.integer(as.character(knn(train, x, cl))))

} # End(donttest)
}


\keyword{methods}
\keyword{spatial}
