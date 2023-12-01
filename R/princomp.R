

setMethod("princomp", signature(x="SpatRaster"),
	function(x, cor=FALSE) {
		stopifnot(nlyr(x) > 1)
		xcov <- layerCor(x, fun="cov", na.rm=TRUE)
		model <- princomp(covmat = xcov$covariance)
		model$center <- diag(xcov$mean)
		if (cor) {
		## Calculate scale as population sd like in in princomp
			S <- diag(xcov$covariance)
			n <- diag(xcov$n)
			model$scale <- sqrt(S * (n-1) / n)
		}
		model
	}
)

