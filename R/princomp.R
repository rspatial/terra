

setMethod("princomp", signature(x="SpatRaster"),
	function(x, cor=FALSE, nobs=FALSE) {
		stopifnot(nlyr(x) > 1)
		xcov <- layerCor(x, fun = "cov", na.rm = TRUE)
		model  <- princomp(covmat = xcov$covariance)
		model$center <- diag(xcov$mean)
		model$n.obs  <- ncell(x) - global(anyNA(x), sum)$sum
		if (cor) {
		## Calculate scale as population sd like in in princomp
			model$n.obs  <- ncell(x) - global(anyNA(x), sum)$sum
			S <- diag(xcov$covariance)
			model$scale <- sqrt(S * (model$n.obs-1) / model$n.obs)
		} else if (nobs) {
			model$n.obs  <- ncell(x) - global(anyNA(x), sum)$sum
		}
		model
	}
)

