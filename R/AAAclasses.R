# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# License GPL v3


setClass("SpatRaster",
	representation (
		ptr = "C++Object"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		if (is.null(object@ptr) || is(object@ptr, "Rcpp_SpatRaster")){
			return(TRUE)
		} else {
			return(FALSE)		
		}
	}
)


setClass("SpatVector",
	representation (
		ptr = "C++Object"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		if (is.null(object@ptr) || is(object@ptr, "Rcpp_SpatVector")){
			return(TRUE)
		} else {
			return(FALSE)		
		}
	}
)


setClass("SpatExtent",
	representation (
		ptr = "C++Object"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		if (is.null(object@ptr) || is(object@ptr, "Rcpp_SpatExtent")){
			return(TRUE)
		} else {
			return(FALSE)		
		}
	}
)


setClass("SpatOptions",
	representation (
		ptr = "C++Object"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		if (is.null(object@ptr) || is(object@ptr, "Rcpp_SpatOptions")){
			return(TRUE)
		} else {
			return(FALSE)		
		}
	}
)

