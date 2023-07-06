# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# License GPL v3


setClass("SpatRaster",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatRaster")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterDataset",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatRasterStack")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterCollection",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatRasterCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVector",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatVectorProxy",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatVectorProxy")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVectorCollection",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatVectorCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatExtent",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatExtent")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatOptions",
	representation (
		pnt = "C++Object"
	),
	prototype (
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatOptions")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatGraticule",
	representation (
		pnt = "C++Object",
		box = "C++Object"
	),
	prototype (
		pnt = NULL,
		pnt = NULL
	),
	validity = function(object)	{
		if (is.null(object@pnt) || is(object@pnt, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)
