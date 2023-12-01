

setMethod("princomp", signature(x="SpatRaster"),
	function(x, cor=FALSE, fix_sign=TRUE) {
		stopifnot(nlyr(x) > 1)
		xcov <- layerCor(x, fun="cov", na.rm=TRUE)
		model <- princomp(covmat = xcov$covariance, cor=cor, fix_sign=fix_sign)
		model$center <- diag(xcov$mean)
		n <- diag(xcov$n)
		if (cor) {
		## Calculate scale as population sd like in in princomp
			S <- diag(xcov$covariance)
			model$scale <- sqrt(S * (n-1) / n)
		}
		model$n <- n
		model
	}
)

