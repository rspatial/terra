

setMethod("princomp", signature(x="SpatRaster"),
	function(x, cor=FALSE) {
		stopifnot(nlyr(x) > 1)
		xcov <- layerCor(x, fun="cov", na.rm=TRUE)
		model <- princomp(covmat = xcov$covariance)
		model$center <- diag(xcov$mean)
		if (cor) {
		## Calculate scale as population sd like in in princomp
			aNA <- anyNA(x)
			x <- mask(x, aNA, maskvalue=TRUE)
			model$n.obs <- ncell(x) - global(aNA, sum)$sum
			n.obs  <- ncell(x) - global(anyNA(x), sum)$sum
			S <- diag(xcov$covariance)
			model$scale <- sqrt(S * (model$n.obs-1) / model$n.obs)
		}
		model
	}
)

