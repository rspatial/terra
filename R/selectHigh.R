
setMethod("selectHighest", signature(x="SpatRaster"),
	function(x, n, low=FALSE) {
		x <- x[[1]]
		n <- min(ncell(x), max(1, n))
		i <- order(values(x, mat=FALSE), decreasing=!low)[1:n]
		x <- rast(x)
		x[i] <- 1
		x
	}
)

