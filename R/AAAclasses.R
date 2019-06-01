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
		return(TRUE)
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
		return(TRUE)
#		object@ptr$valid
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
		object@ptr$valid
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
		return (true)
	}
)

