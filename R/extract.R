
setMethod("extract", signature(x="SpatRaster", y="SpatLayer"), 
function(x, y, fun="", ...) { 
    r <- x@ptr$extractLayer(y@ptr, fun)
	x <- show_messages(x, "extract")		
	r
}
)


setMethod("[", c("SpatRaster", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	theCall <- sys.call(-1)
	narg <- length(theCall) - length(match.call(call=sys.call(-1)))
	if (narg > 0) {
		i <- cellFromRow(x, i)
	} 
	r <- x@ptr$extractCell(i)
	show_messages(x)
	r
})


