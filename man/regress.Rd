\name{regress}

\docType{methods}

\alias{regress}
\alias{regress,SpatRaster,SpatRaster-method}
\alias{regress,SpatRaster,numeric-method}
\alias{regress,SpatRaster,data.frame-method}

\title{Cell level regression}

\description{
Run a regression model for each cell of a SpatRaster. The independent variable can be defined either by:
\itemize{
  \item a \code{numeric} vector of length \code{nlyr(y)} -- one covariate that varies across the layers but is constant across cells (e.g. time);
  \item a \code{SpatRaster} of the same number of layers as \code{y} -- one covariate that varies both across layers and across cells;
  \item a \code{data.frame} with one row per layer of \code{y} and one column per predictor -- multiple covariates that vary across layers but are constant across cells (e.g. an experimental design like NPK fertilizer levels).
}
}

\usage{

\S4method{regress}{SpatRaster,numeric}(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

\S4method{regress}{SpatRaster,SpatRaster}(y, x, formula=y~x, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)

\S4method{regress}{SpatRaster,data.frame}(y, x, formula=NULL, na.rm=FALSE, cores=1, filename="", overwrite=FALSE, ...)
}


\arguments{
\item{y}{SpatRaster. The dependent variable. Each layer is one observation per cell.}
\item{x}{The independent variable(s):
\itemize{
  \item \code{SpatRaster} of the same number of layers as \code{y};
  \item \code{numeric} vector of length \code{nlyr(y)};
  \item \code{data.frame} with \code{nlyr(y)} rows and one column per predictor (no \code{NA}s allowed).
}}
\item{formula}{regression formula. For the \code{numeric} and \code{SpatRaster} methods the default is \code{y ~ x}; you can add additional terms such as \code{I(x^2)}. For the \code{data.frame} method the default (when \code{formula=NULL}) is \code{y ~ .}, expanding to all columns of \code{x}; you can also pass a custom formula whose right-hand side references the column names of \code{x} (the left-hand side is ignored).}
\item{na.rm}{logical. If \code{TRUE}, layers with \code{NA} response in a given cell are dropped before fitting. A cell whose remaining row count drops below the number of regression coefficients yields all \code{NA}s.}
\item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object.}
\item{filename}{character. Output filename}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster. One layer per regression coefficient, named after the corresponding column of \code{model.matrix(formula, x)} (e.g. \code{(Intercept)}, the names of the predictors, and any interaction or polynomial terms).}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))
x <- regress(s, 1:nlyr(s))

# data.frame method: a small NPK fertilizer experiment.
# Suppose `y` has one layer per (N, P, K) treatment.
\dontrun{
NPK <- data.frame(
    N = c(0,  50, 100,   0,  50, 100,   0,  50, 100),
    P = c(0,   0,   0,  25,  25,  25,  50,  50,  50),
    K = c(0,  10,  20,   0,  10,  20,   0,  10,  20)
)
# Default formula y ~ N + P + K:
fit <- regress(y, NPK)

# Custom formula with interactions and a quadratic term:
fit2 <- regress(y, NPK, formula = yield ~ N + P + K + I(N^2) + N:P)
}
}

\keyword{methods}
\keyword{spatial}

