# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# License GPL v3


setClass("SpatRaster",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatRaster")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterDataset",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatRasterStack")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatRasterCollection",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatRasterCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVector",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)

setClass("SpatVectorProxy",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatVectorProxy")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatVectorCollection",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatVectorCollection")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatExtent",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatExtent")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)


setClass("SpatOptions",
	representation (
		pntr = "C++Object"
	),
	prototype (
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatOptions")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)



setClass("SpatGraticule",
	representation (
		pntr = "C++Object",
		box = "C++Object"
	),
	prototype (
		pntr = NULL,
		pntr = NULL
	),
	validity = function(object)	{
		if (is.null(object@pntr) || is(object@pntr, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
)
