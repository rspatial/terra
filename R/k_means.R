
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



h_clust <- function(x, ngroups, dist_metric="euclidean", clust_method="complete", agfun=mean, matchfun="squared", ..., maxcell=10000, filename="", overwrite=FALSE, wopt=list()) {
	stopifnot(maxcell > 0)
	stopifnot(ngroups > 0)
	stopifnot(ngroups < maxcell)
	d <- na.omit(spatSample(x, maxcell, "regular"))
	dd <- stats::dist(d, distmetric)
	hc <- stats::hclust(dd, clustmethod)
	th <- sort(hc$height, TRUE)[ngroups]
	cls <- stats::cutree(hc, h = th)
	hc <- cut(stats::as.dendrogram(hc), h=th)$upper
	d <- aggregate(d, list(cls=cls), agfun)
	cls <- d$cls 
	d$cls <- NULL
	b <- bestMatch(x, d, fun=matchfun, ..., filename=filename, overwrite=overwrite, wopt=wopt)	
	return(list(clusters=b, dendrogram=hc))
}

