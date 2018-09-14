# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# Licence GPL v3



setClass('SpatPolygons',
	representation (
		ptr = "Rcpp_SpatPolygons"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
#		object@ptr$valid
	}
)


setClass('SpatExtent',
	representation (
		ptr = "Rcpp_SpatExtent"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		object@ptr$valid
	}
)



setClass('SpatRaster',
	representation (
		ptr = "Rcpp_SpatRaster"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		return(TRUE)
	}
)



