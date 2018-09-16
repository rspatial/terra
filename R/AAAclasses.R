# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : November 2017
# Version 0
# Licence GPL v3



setClass('SpatPolygons',
	representation (
		ptr = "C++Object"
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
		ptr = "C++Object"
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
		ptr = "C++Object"
	),	
	prototype (	
		ptr = NULL
	),
	validity = function(object)	{
		return(TRUE)
	}
)



