\name{predict}

\docType{methods}

\alias{predict}
\alias{predict,SpatRaster-method}

\title{Spatial model predictions}

\description{
Make a SpatRaster object with predictions from a fitted model object (for example, obtained with \code{lm}, \code{glm}). The first argument is a SpatRaster object with the predictor variables. The \code{\link{names}} in the Raster object should exactly match those expected by the model. Any type of model (e.g. glm, gam, randomForest) for which a predict method has been implemented (or can be implemented) can be used. 

This approach of using model predictions is commonly used in remote sensing (for the classification of satellite images) and in ecology, for species distribution modeling.
}

\usage{
\S4method{predict}{SpatRaster}(object, model, fun=predict, ..., filename="", format="", datatype="FLT4S", overwrite=FALSE)
}

\arguments{
  \item{object}{SpatRaster}
  \item{model}{fitted model of any class that has a }predict" method (or for which you can supply a similar method as \code{fun} argument. E.g. glm, gam, or randomForest}
  \item{fun}{function. Default value is \code{predict}, but can be replaced with e.g. predict.se (depending on the type of model), or your own custom function}
  \item{...}{Additional arguments for \code{fun}}    

  \item{filename}{Character. Output filename. Optional}
  \item{format}{Character. Output file format. See \code{\link{writeRaster}}}
  \item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
 }


\value{
SpatRaster
}

\keyword{methods}
\keyword{spatial}
