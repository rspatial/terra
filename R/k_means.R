
setMethod("k_means", signature(x="ANY"),
	function(x, centers=3, ...) {
		stats::kmeans(x, centers=centers, ...)
	}
)

setMethod("k_means", signature(x="SpatRaster"),
	function(x, centers=3, ..., maxcell=1000000, filename="", overwrite=FALSE, wopt=list()) {
		stopifnot(maxcell > 0)
		if (ncell(x) <= maxcell) {
			v <- na.omit(as.matrix(x))
			omit <- as.vector(attr(v, "na.action"))
			km <- kmeans(v, centers=centers, ...)
			r <- rast(x, nlyr=1)
			if (is.null(omit)) {
				values(r) <- km$cluster
			} else {
				r[-omit] <- km$cluster
			}
			if (filename != "") {
				r <- writeRaster(r, filename=filename, overwrite=overwrite, wopt=wopt)
			}
		} else {
			pkmeans = function(x, newdata) {
				apply(newdata, 1, function(i) which.min(colSums((t(x$centers) - i)^2)))
			}
			v <- na.omit(spatSample(x, maxcell, "regular"))
			km <- kmeans(v, centers=centers, ...)
			r <- predict(logo, km, fun=pkmeans, na.rm=TRUE, filename=filename, overwrite=overwrite, wopt=wopt)
		}
		r
	}
)

