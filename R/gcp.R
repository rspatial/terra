

setClass("GCP",
	representation (
		gcp = "matrix"
	),
	prototype (
		gcp = cbind(fx=0, fy=0, tx=0, ty=0)[0,]
	)
)

if (!isGeneric("addGCP")) {setGeneric("addGCP", function(x, ...) standardGeneric("addGCP"))}

setMethod("show", signature(object="GCP"), 
	function(object) {
		m <- object@gcp
		show(m)
		if (!is.null(dev.list())) {
			for (i in 1:nrow(m)) {
				arrows(m[i,1], m[i,2], x1 = m[i,3], y1 = m[i,4], col="red", length = 0.1)
			}
		}
	}
)


setMethod("addGCP", signature(x="GCP"), 
	function(x, from_to) {
		if (missing(from_to)) {
			from_to = terra:::RS_locator(2, "l")
			from_to <- rbind(as.vector(t(from_to)))
		} 
		if (ncol(from_to) == 4) {
			x@gcp <- rbind(x@gcp, from_to)
		}
		arrows(from_to[1,1], from_to[1,2], x1 = from_to[1,3], y1 = from_to[1,4], col="red", length = 0.1)
		x
	}
)


#gcp <- new("GCP")
#gcp <- addGCP(gcp)
