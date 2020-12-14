\name{interpolate}

\docType{methods}

\alias{interpolate}
\alias{interpolate,SpatRaster-method}

\title{Interpolate}

\description{
Make a SpatRaster with interpolated values using a fitted model object of classes such as "gstat" (gstat package) or "Krige" (fields package). That is, these are models that have location ("x" and "y", or "longitude" and "latitude") as independent variables. If x and y are the only independent variables provide an empty (no associated data in memory or on file) SpatRaster for which you want predictions. If there are more spatial predictor variables provide these as a SpatRaster in the first argument of the function. If you do not have x and y locations as implicit predictors in your model you should use \code{\link[terra]{predict}} instead.
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
\item{const}{data.frame. Can be used to add a constant for which there is no SpatRaster for model predictions. This is particularly useful if the constant is a character-like factor value}
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


\keyword{methods}
\keyword{spatial}
