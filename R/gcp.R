

setClass("GCP",
	representation (
		gcp = "matrix"
	),
	prototype (
		gcp = cbind(fx=0, fy=0, tx=0, ty=0)[0,]
	)
)



setMethod("show", signature(object="GCP"),
	function(object) {
		m <- object@gcp
		show(m)
		if (!is.null(grDevices::dev.list())) {
			for (i in 1:nrow(m)) {
				graphics::arrows(m[i,1], m[i,2], x1 = m[i,3], y1 = m[i,4], col="red", length = 0.1)
			}
		}
	}
)


#setMethod("add<-", signature(x="GCP"),
#	function(x, value) {
#		if (missing(value)) {
#			value <- terra:::RS_locator(2, "l")
#			value <- rbind(as.vector(t(value)))
#		}
#		if (ncol(value) == 4) {
#			x@gcp <- rbind(x@gcp, value)
#		}
#		if (!is.null(grDevices::dev.list())) {
#			graphics::arrows(value[1,1], value[1,2], x1 = value[1,3], y1 = value[1,4], col="red", length = 0.1)
#		}
#		x
#	}
#)


#gcp <- new("GCP")
#gcp <- addGCP(gcp)
