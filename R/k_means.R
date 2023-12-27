
setMethod("k_means", signature(x="ANY"),
	function(x, centers=3, ...) {
		stats::kmeans(x, centers=centers, ...)
	}
)

setMethod("k_means", signature(x="SpatRaster"),
	function(x, centers=3, ..., maxcell=1000000, filename="", overwrite=FALSE, wopt=list()) {
		stopifnot(maxcell > 0)
		if (ncell(x) <= maxcell) {
			v <- na.omit(values(x))
			omit <- as.vector(attr(v, "na.action"))
			km <- stats::kmeans(v, centers=centers, ...)
			out <- rast(x, nlyr=1)
			if (is.null(omit)) {
				values(out) <- km$cluster
			} else {
				out[-omit] <- km$cluster
			}
			if (filename != "") {
				out <- writeRaster(out, filename=filename, overwrite=overwrite, wopt=wopt)
			}
		} else {
			#pkmeans = function(x, newdata) {
			#	apply(newdata, 1, function(i) which.min(colSums((t(x$centers) - i)^2)))
			#}

			pkmeans <- function(x, newdata) {
				vec <- integer(nrow(newdata))
				newdata <- as.matrix(newdata)
				for (i in seq_len(nrow(newdata))) {
					vec[i] <- which.min(colSums((t(x) - newdata[i, ])^2))
				}
				vec
			}			
			v <- unique(na.omit(spatSample(x, maxcell, "regular")))
			km <- stats::kmeans(v, centers=centers, ...)$centers		
			out <- predict(x, km, fun=pkmeans, na.rm=TRUE, filename=filename, overwrite=overwrite, wopt=wopt)
		}
		out
	}
)

