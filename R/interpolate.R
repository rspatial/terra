# Author: Robert J. Hijmans
# Date :  August 2009
# Version 0.9
# License GPL v3


setMethod("interpolate", signature(object="SpatRaster"),
	function(object, model, fun=predict, ..., xyNames=c("x", "y"), factors=NULL, const=NULL, index=NULL, cores=1, cpkgs=NULL, na.rm=FALSE, filename="", overwrite=FALSE, wopt=list()) {

		out <- rast(object)
		hv <- hasValues(object)
		nms <- c(xyNames, names(object))
		if (length(unique(nms)) != length(nms)) {
			tab <- table(nms)
			error("interpolate", "duplicate names: ", tab[tab>1])
		}
		nc <- ncol(out)
		testrow <- round(0.51*nrow(object))
		ntest <- min(nc, 500)
		xy <- xyFromCell(out, cellFromRowCol(out, testrow, 1):cellFromRowCol(out, testrow, ntest))
		colnames(xy) <- xyNames
		if (hv) {
			readStart(object)
			on.exit(readStop(object))
			d <- readValues(object, testrow, 1, 1, ntest, TRUE, TRUE)
			xy <- cbind(xy, d)
		}
		r <- .runModel(model, fun, xy, 1, const, (na.rm & hv), index, cores=NULL, ...)
		rdim <- dim(r)
		if (!is.null(rdim)) {
			if (rdim[1] == 1) {
				nl <- rdim[1]
			} else {
				nl <- rdim[2]
			}
		} else {
			nl <- 1
		}
		out <- rast(object, nlyrs=nl)
		cn <- colnames(r)
		if (length(cn) == nl) names(out) <- make.names(cn, TRUE)

		doclust <- FALSE
		if (inherits(cores, "cluster")) {
			doclust <- TRUE
		} else if (cores > 1) {
			doclust <- TRUE
			cores <- parallel::makeCluster(cores)
			on.exit(parallel::stopCluster(cores), add=TRUE)
		}
		if (doclust) {
			parallel::clusterExport(cores, c("model", "fun"), environment())
			if (!is.null(cpkgs)) {
				parallel::clusterExport(cores, "cpkgs", environment())
				parallel::clusterCall(cores, function() for (i in 1:length(cpkgs)) {library(cpkgs[i], character.only=TRUE) })
			}
			export_args(cores, ..., caller="interpolate")			
		} else {
			cores <- NULL
		}
		
		b <- writeStart(out, filename, overwrite, sources=sources(object), wopt=wopt)
		for (i in 1:b$n) {
			xy <- xyFromCell(out, cellFromRowCol(out, b$row[i], 1):cellFromRowCol(out, b$row[i]+b$nrows[i]-1, nc))
			colnames(xy) <- xyNames
			if (hv) {
				d <- readValues(object, b$row[i], b$nrows[i], 1, nc, TRUE, TRUE)
				xy <- cbind(xy, d)
			}
			v <- .runModel(model, fun, xy, nl, const, (na.rm & hv), index, cores=cores, ...)
			writeValues(out, v, b$row[i], b$nrows[i])
		}
		writeStop(out)
	}
)


