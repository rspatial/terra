
setMethod("extract", signature(x="SpatRaster", y="SpatLayer"), 
function(x, y, fun="", ...) { 
    r <- x@ptr$extractLayer(y@ptr, fun)
	x <- show_messages(x, "extract")		
	r
})


setMethod("[", c("SpatRaster", "missing", "missing"),
function(x, i, j, ... , drop=TRUE) {
	values(x)
})


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
#	if (any(i) > 2^.Machine$double.digits) { warning() }
	if (nargs() == 3) {
		i <- cellFromRowCol(x, i, 1:ncol(x))
	} 
	r <- x@ptr$extractCell(i)
	show_messages(x)
	r
})


setMethod("[", c("SpatRaster", "numeric", "numeric"),
function(x, i, j, ... ,drop=TRUE) {
	i <- cellFromRowCol(x, i, j)
#	if (any(i) > 2^.Machine$double.digits) { warning() }
	r <- x@ptr$extractCell(i)
	show_messages(x)
	r
})

