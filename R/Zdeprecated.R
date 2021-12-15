
## spatstat conflicts

if (!isGeneric("area")) {setGeneric("area", function(x, ...) standardGeneric("area"))}

setMethod("area", signature(x="SpatRaster"), 
	function(x, sum=TRUE, correct=FALSE, mask=FALSE, filename="", ...) {
		if (!sum) {
			error("area", 'area(x, sum=FALSE) will be removed. Use "cellSize(x)"')
		} else {
			error("area", 'area(x, sum=TRUE) will be removed. Use "expanse(x)" or "global(cellSize(x), "sum")"')
		}
	}
)

setMethod("area", signature(x="SpatVector"), 
	function(x, ...) {
		error("area", 'area was removed. Use "expanse(x)"')
	}
)


