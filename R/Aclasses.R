# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# License GPL v3


setClass("SpatRaster",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatRaster")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterDataset",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatRasterStack")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterCollection",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatRasterCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVector",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatVectorProxy",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatVectorProxy")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVectorCollection",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatVectorCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatExtent",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatExtent")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatOptions",
	representation (
		cpp = "C++Object"
	),
	prototype (
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatOptions")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatGraticule",
	representation (
		cpp = "C++Object",
		box = "C++Object"
	),
	prototype (
		cpp = NULL,
		cpp = NULL
	),
	validity = function(object)	{
		if (is.null(object@cpp) || is(object@cpp, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)
