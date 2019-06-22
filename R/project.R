
if (!isGeneric("project")) {setGeneric("project", function(x,...) standardGeneric("project"))}	

setMethod("project", signature(x="SpatVector"), 
	function(x, crs, ...)  {
		x@ptr <- x@ptr$transform_crs(crs)
		show_messages(x, "project")
	}
)
